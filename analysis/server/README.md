# Running the models on Artemis (Sussex HPC)

Model fitting is too heavy for a laptop, so it runs as SLURM array jobs on
**Artemis**. The `./hpc` script in this folder wraps the whole loop — push
code, submit, watch, pull results — over SSH.

**One job fits one model.** `models.R` is the registry of what a model *is*;
`./hpc fit <model>` submits an array for it, each array task writing one shard;
`./hpc combine <model>` merges that model's shards into a single fit.

## The cluster, in a browser

Everything below drives Artemis over SSH, but the web interface is useful for
looking at files, checking a job by hand, and installing your SSH key the first
time. All of it needs the **GlobalProtect VPN** connected first.

| | |
| --- | --- |
| Open OnDemand (OOD) | <https://ood.artemis.hrc.sussex.ac.uk/> |
| a shell on the login node | OOD → Clusters → `>_ artemis Shell Access` |
| submit a `.slurm` by hand | OOD → Jobs → Jobs Composer |
| browse your files | `https://ood.artemis.hrc.sussex.ac.uk/pun/sys/dashboard/files/fs//mnt/lustre/users/psych/<user>/IGComputational` |
| Artemis documentation | <https://artemis-docs.hpc.sussex.ac.uk/artemis/> |

Storage, for orientation:

| path | use |
| --- | --- |
| `/mnt/lustre/users/<group>/<user>/` | long-term — code and fitted models live here |
| `/mnt/lustre/scratch/<group>/<user>/` | fast scratch — job logs go here |
| `/mnt/nfs2/<group>/<user>/` | home directory |

## One-time setup (per machine)

```bash
bash analysis/server/setup-ssh.sh
```

Generates a machine-local keypair at `~/.ssh/artemis`, adds the `artemis`
alias to `~/.ssh/config`, and prints the exact command to install the public
key on the cluster (OOD -> Clusters -> `>_ artemis Shell Access`). Idempotent.

Run it on every machine you work from — each gets its own key, and Artemis
accepts as many as you append to `~/.ssh/authorized_keys`. Never copy a
private key between machines.

## Every session

The **GlobalProtect VPN must be connected** before anything below works.

```bash
cd analysis/server
./hpc check                # confirms VPN + SSH are good
./hpc push                 # copy *.R and *.slurm to the cluster (CRLF -> LF)
./hpc models               # what can be fitted
./hpc fit gam_lnr          # submit one array job for one model
./hpc queue                # what's running
./hpc log gam_lnr          # tail that model's newest .out/.err
./hpc ls                   # fitted .rds on the cluster
./hpc combine gam_lnr      # merge its shards, add waic, delete the shards
./hpc pull                 # bring combined fits into analysis/models/
```

`push` **mirrors** rather than merges: a top-level `*.R` / `*.slurm` on the
cluster that no longer exists locally is deleted, so a renamed script can never
be submitted by accident. `models/`, `tests/` and the logs are untouched.

Extra arguments to `fit` and `combine` are passed to `sbatch`, so a smoke test
is:

```bash
IGC_MODELS_DIR=/mnt/lustre/users/psych/dmm56/IGComputational/smoke \
IGC_NPARTICIPANTS=30 IGC_WARMUP=300 IGC_SAMPLES=100 \
  ./hpc fit gam_lnr --array=1-2 --partition=short --time=01:00:00 --cpus-per-task=8 --mem=16G
```

Always give a test run its **own `IGC_MODELS_DIR`**. Shards are named
`<model>_<illusion>_<shard>.rds` and `file_refit = "never"` means a shard that
already exists is silently kept — a 30-participant test shard left in the
production directory would be adopted by the production run.

## The production runs

Settings below are the ones measured in `AGENT.md` §4. Full data, cold starts,
`long` partition, one job per model:

```bash
./hpc push
./hpc install cogmod      # only if cogmod moved
./hpc precompile          # only after a toolchain change
./hpc fit gam_lnr
./hpc fit gam_lnr6
./hpc fit gam_ddm4
```

The defaults in `fit.slurm` are that production configuration: `--array=1-4`,
`--cpus-per-task=16`, `--mem=32G`, `--partition=long`, `--time=2-00:00:00`,
with `IGC_NPARTICIPANTS=all`, `IGC_WARMUP=1000`, `IGC_SAMPLES=500` and
`IGC_CHAINS=2`. That is 4 shards x 2 chains x 500 = **4,000 draws per model**.

Three models is 12 tasks against `long`'s 140-CPU per-user cap, which allows
`140 / 16 = 8` at a time, so the last four start as earlier ones finish. That
costs nothing: `--time` is per task, not per array, and pending time does not
count against it.

When they are done:

```bash
./hpc combine gam_lnr
./hpc combine gam_lnr6
./hpc combine gam_ddm4
./hpc pull                # combined/*.rds -> analysis/models/
```

`combine` is one job, not an array, and its constraint is memory rather than
CPU: `--mem=128G` on `general` for 6 h, because `add_criterion("waic")` builds
a 1500 x 323,981 pointwise log-likelihood matrix (~3.9 GB) through an R-level
`log_lik` called once per response.

It strips each shard's `$file` slot before merging. A shard is written by
`brm(file = ...)` and so remembers its own path; `combine_models()` keeps the
first shard's, and `add_criterion()` writes the fit back to whatever `$file`
says — which silently saved the *combined* fit over shard 1, after which the
mini merged that with shard 2 and came out larger than the full fit. Found on
2026-09-18 by running `combine` against real shards for the first time.

## Models

`models.R` holds one entry per model: the illusion it is fitted to and a
function returning the brms formula. It is the only place a model is defined —
`fit_model.R` fits the one named by `IGC_MODEL`, `combine_model.R` merges the
same one, and `./hpc` reads the *names* straight out of the file (which is why
the declaration lines must stay in the form `  <name> = list(`).

| model | family | distributional parameters with a 2-D smooth |
| --- | --- | --- |
| `gam_lnr` | `cogmod_lnr()` | drift, `nuone`, `sigmazero`, `sigmaone`, `ndt` (`sigmabias = 0`) |
| `gam_lnr6` | `cogmod_lnr()` | as LNR plus `sigmabias`, the between-trial start-point range |
| `gam_ddm4` | `cogmod_ddm()` | drift, `boundary`, `bias`, `ndt` (all three between-trial SDs 0) |
| `gam_ddm5` | `cogmod_ddm()` | as DDM-4 plus `sigmadrift` |

Every smooth is `t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs =
c("cr", "cr")) + (1 | Participant)`; `poutlier` is `1 + (1 | Participant)`.
Production is `gam_lnr`, `gam_lnr6` and `gam_ddm4`; `gam_ddm5` is defined and
submittable but is not part of the final run. Note that `gam_lnr6` should cost
*more* per gradient than `gam_lnr`, not less: a free `sigmabias` takes cogmod's
erfc-based two-tail path instead of the single-tail shortcut, which 0.3.3
measured at ~15% dearer (and ~20% cheaper for `sigmabias = 0`).

Adding a model is one entry in `models.R` and nothing else.

## Paths

Nothing is hard-coded to one account. `./hpc` passes the paths to `sbatch`
via `--output/--error/--chdir/--export`, and the scripts read them from
`IGC_USERS_DIR` / `IGC_MODELS_DIR` / `IGC_SCRATCH_DIR`. Defaults:

| variable | default |
| --- | --- |
| `IGC_HPC_USER` | `dmm56` |
| `IGC_HPC_GROUP` | `psych` |
| `IGC_USERS_DIR` | `/mnt/lustre/users/$IGC_HPC_GROUP/$IGC_HPC_USER/IGComputational` |
| `IGC_SCRATCH_DIR` | `/mnt/lustre/scratch/$IGC_HPC_GROUP/$IGC_HPC_USER/IGComputational` |
| `IGC_MODELS_DIR` | `$IGC_USERS_DIR/models` |

`IGC_HPC_USER` and `IGC_HPC_GROUP` rewrite all four paths at once, which is
what makes a second account cheap — see the next section. Note it changes the
*paths*, not who you log in as: that comes from the `artemis` entry in your own
`~/.ssh/config`. Pointing at someone else's directories while logged in as
yourself gets you permission denied, not their run.

Submitting a `.slurm` file by hand (e.g. from the OOD Jobs Composer) still
works — the scripts then fall back to `$SLURM_SUBMIT_DIR` — but you must set
`IGC_MODEL` yourself, since that is what selects the model.

## Running from a second cluster account

The per-user quota is the binding constraint (`long` allows one user 140 CPUs,
i.e. 8 of our tasks at a time), and it is **per account**. A colleague with
their own Artemis account and a clone of this repo can fit a different subset
of the models at the same time, roughly doubling throughput. The registry
travels in git; only the account is local.

Their setup, once:

```bash
IGC_HPC_USER=oc236 bash analysis/server/setup-ssh.sh   # their key, their ssh alias
echo 'IGC_HPC_USER=oc236' > analysis/server/hpc.local   # gitignored; no tracked file changes
cd analysis/server
./hpc check
./hpc setup          # create their project dirs
./hpc install        # build their R library -- and CmdStan, ~25 min the first time
./hpc precompile     # their own CmdStan precompiled header
./hpc push
```

`hpc.local` is sourced by `./hpc` before any default is applied and is in
`.gitignore`, so it never reaches a commit. Add `IGC_HPC_GROUP=informatics`
(or whatever) if they are not in the `psych` tree, or set `IGC_USERS_DIR` /
`IGC_SCRATCH_DIR` outright if their layout differs. Everything else — module,
library path, partitions, the `IGC_*` run-shaping variables — is identical.

Then split the models by name. Both sides have the same `models.R`, so the only
thing to agree on is who runs what:

```bash
# dmm56
./hpc fit gam_lnr ; ./hpc fit gam_lnr6 ; ./hpc fit gam_ddm4
# oc236
./hpc fit gam_ddm5 ; ./hpc fit <whatever is added next>
```

A model that only one of you fits still has to be **defined in the shared
`models.R` and committed**, or `./hpc fit` will reject the name. That is the
point: the definition is reviewed and version-controlled, the account is not.

### Collecting the results

Lustre home directories are world-readable within the cluster (`drwxr-xr-x`,
and the fits themselves `-rw-r--r--`), so whoever is assembling the analysis
can pull the other account's combined fits directly, without either of them
copying anything by hand:

```bash
# from dmm56's machine, after oc236 has run ./hpc combine
IGC_MODELS_DIR=/mnt/lustre/users/psych/oc236/IGComputational/models ./hpc pull
```

This reads over SSH as *your* account and writes into your local
`analysis/models/`, which is gitignored. Ask them to run `./hpc combine
<model>` first — a directory of raw shards is not what you want, and the
shards are deleted once combined.

If they would rather hand the files over than have you read their directory,
`./hpc pull` on their side puts the same files in their own
`analysis/models/`.

## Run-shaping variables

`./hpc fit` forwards these to the job when you set them; everything else keeps
the script default.

| variable | default | purpose |
| --- | --- | --- |
| `IGC_NPARTICIPANTS` | `all` | subset size for tests, e.g. `30` |
| `IGC_WARMUP` | `1000` | warmup iterations |
| `IGC_SAMPLES` | `500` | post-warmup draws per chain |
| `IGC_CHAINS` | `2` | chains per array task; threads per chain is `cpus / chains` |
| `IGC_FILE_REFIT` | `never` | `always` forces a clean refit |
| `IGC_KEEP_SHARDS` | unset | set to keep shards after combining |
| `IGC_COGMOD_REF` | `dev` | which cogmod branch/tag `./hpc install` tracks |

Each fit prints a `REPORT ...` line (wall time, mean leapfrog steps,
treedepth, step size, divergences, max Rhat, min ESS ratio) to its `.out`, so
runs can be compared without pulling the `.rds`.

## Precompile the CmdStan header before a cold array

Run this once after any change to the toolchain, the cmdstan version, or the
`stan_model_args` in `fit_model.R`:

```bash
./hpc precompile     # then
./hpc fit gam_lnr
```

CmdStan precompiles a Stan header into
`~/.cmdstan/<version>/stan/src/stan/model/model_header.hpp.gch/`. Note that
this is a **directory**, not a file: it holds one ~800 MB variant per
compiler/flag combination, e.g.

```
model_header_12_3.hpp.gch                    # GCC 12.3, plain
model_header_nochecks_12_3.hpp.gch           # GCC 12.3, no range checks
model_header_threads_nochecks_12_3.hpp.gch   # GCC 12.3, threaded, no range checks
```

Our jobs need the last of those (foss-2023a = GCC 12.3, `threading()`,
`STAN_NO_RANGE_CHECKS`) — and note that a run *without* threading builds a
different variant, so having *a* `.gch` there is not the same as having ours.
If the variant we need is missing when an array starts, **every task races to
build the same file**, and all but one die with:

```
stan/src/stan/model/model_header.hpp:2:39:
    error: while reading precompiled header: No such file or directory
Error: An error occured during compilation!
```

This killed 10 of 12 tasks in job 11366858 on 2026-09-17, in ~3.5 minutes
each. Worse, the loser tasks leave a **corrupt** variant behind, after which
nothing compiles at all — not even a trivial model — until CmdStan is
rebuilt:

```bash
./hpc sh "... Rscript -e 'cmdstanr::rebuild_cmdstan(cores = 8)'"
```

Two things that do *not* work as a shortcut:

- Adding `PRECOMPILED_HEADERS=false` to `~/.cmdstan/<version>/make/local` —
  cmdstanr rewrites `make/local` from its own `cpp_options` on every compile,
  so the setting is discarded.
- Deleting `stan/src/stan/model/*.gch` — that glob matches the *directory*,
  and `rm -f` refuses it ("Is a directory"). The stale variants survive.

So: keep `precompile.R`'s `cpp_options` identical to `stan_model_args` in
`fit_model.R`, and run `./hpc precompile` before a cold array.

## Partitions and resource limits

Per-user quotas, confirmed live with
`sacctmgr show qos format=Name,MaxTRESPU,MaxWall`:

| partition | max runtime | max CPUs | max RAM |
| --- | --- | --- | --- |
| `short` | 2 hours | 550 | 3.6 TB |
| `general` (default) | 8 hours | 400 | 2.7 TB |
| `long` | 8 days | **140** | 900 GB |
| `verylong` | 30 days | 70 | 900 GB |
| `gpu` | 3 days | 300 | 2.1 TB |

There is **no 24-hour tier**, and the 3-day `gpu` partition is GPU-only by
policy. Shorter runtime buys a bigger allowance. Defaults if you ask for
nothing are 1 CPU per task, 4 GB RAM per CPU, and the `general` partition.

**This directly caps concurrency**, at `floor(140 / cpus-per-task)` on `long`.
With `--cpus-per-task=16` that is 8 tasks at once across *all* your jobs —
which is why production is 4 tasks per model for two models rather than 8 for
one. Anything beyond it sits in `PENDING (QOSMaxCpuPerUserLimit)`.

Memory is not the binding constraint: 8 x 32 GB = 256 GB against `long`'s
900 GB ceiling.

### Getting jobs dispatched sooner

Slurm priority is dominated by **Partition** and **Association** (both "high"),
then Age and TRES. JobSize is inversely proportional to the request, so
slimmer jobs start sooner. Practical consequences:

- **Always set `--time` explicitly** rather than inheriting the partition
  default. Reserving 8 days for a 20-hour job inflates predicted consumption
  and delays dispatch. `fit.slurm` sets `2-00:00:00`; override with
  `./hpc fit gam_lnr --time=12:00:00`.
- Request only the CPUs/RAM actually used — over-requesting blocks resources
  and enlarges your apparent job size.
- Many small tasks beat one huge one.
- Avoid sustained bursts of high-volume processing; FairShare penalises it
  (idle interactive sessions hurt most).

Check where you stand with `sprio -u dmm56`, or:

```bash
./hpc sh "squeue --Format=JobID,State,Reason,PriorityLong,Partition -u dmm56"
```

Jobs longer than a partition's limit are killed with a message in the log. Stan
sampling cannot checkpoint mid-chain — our equivalent is the array itself: each
task writes its own `.rds`, so a lost task costs one shard rather than the
whole run, and resubmitting skips the shards that finished.

## R environment on the cluster

The jobs load **one module** that already provides most of the stack:

```
CmdStanR/0.7.1-foss-2023a-R-4.3.2   # R 4.3.2 + mgcv + brms 2.21.0 + cmdstanr 0.7.1 + dplyr
```

This is the only module combination on Artemis that ships `mgcv`, `brms` and
`cmdstanr` together. `R/4.4.1-gfbf-2023b` loads, but has no `mgcv` and no
`cmdstanr`, so the `t2()` smooths and the cmdstanr backend both fail there.

What the module does *not* ship (`datawizard`, `cogmod`, and a newer `cmdstanr`)
lives in a project library, built by `install_pkgs.R` via `./hpc install`:

```
/mnt/lustre/users/psych/$IGC_HPC_USER/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3
```

### Refreshing cogmod

`cogmod` is installed from GitHub, so it needs an explicit reinstall to pick up
new commits. Naming a package on the command line force-reinstalls it:

```bash
./hpc install cogmod     # pull cogmod up to the latest commit on IGC_COGMOD_REF
./hpc install            # only install what is missing (or below the floor)
./hpc install all        # force-reinstall everything
```

This runs through `srun` on a compute node rather than compiling on the login
node, and prints the resulting versions plus the cogmod ref and commit SHA at
the end.

**The fits require cogmod >= 0.3.3**, which fixed the non-finite tail gradient
in the LNR and the init jitter that started smooths far from their targets (see
`cogmod_inits_issue.md`). `install_pkgs.R` fails if the installed version is
below that floor, and `fit_model.R` refuses to start — better a failed install
than a 20-hour job with the old numerics. `IGC_COGMOD_REF` is `dev` until 0.3.3
is merged into `main`; change the default in `hpc` when it lands.

CmdStan itself is already built at `~/.cmdstan/cmdstan-2.39.0` on `dmm56` and
is found automatically by `cmdstanr`. On an account that has none, `./hpc
install` builds one (`install_cmdstan()`, ~15-25 min, once per account) before
going on to the packages.

Both the module name and the library path are overridable:

| variable | default |
| --- | --- |
| `IGC_R_MODULE` | `CmdStanR/0.7.1-foss-2023a-R-4.3.2` |
| `IGC_R_LIBS` | `/mnt/lustre/users/psych/$IGC_HPC_USER/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3` |

`./hpc fit` / `./hpc combine` pass both to the job via `--export`; the `.slurm`
scripts fall back to these same defaults when submitted by hand.

Two Artemis quirks the `.slurm` files work around:

1. SLURM runs batch scripts in a *non-interactive* shell, where `module` is
   not defined — so they source `/etc/profile.d/lmod.sh` first.
2. Compute nodes get only `/opt/ohpc/pub/modulefiles` on `MODULEPATH`; the
   EasyBuild tree holding R and CmdStanR is on the login node's path only.
   The scripts therefore `module use /mnt/shared/easybuild/modules/all`
   (overridable with `IGC_EB_MODULES`) before `module load`. Without this a
   job dies with "module(s) exist but cannot be loaded as requested".

Data is read straight from the GitHub raw URLs in `fit_model.R` — the compute
nodes do have outbound internet, so nothing needs pushing but code.

`.gitattributes` forces LF on `*.slurm` so a CRLF file never reaches the
cluster (which fails with `bad interpreter`). `./hpc push` strips CR as well,
belt and braces.

## If SSH starts refusing connections

`kex_exchange_identification: read: Connection reset` means sshd is
rate-limiting, not that the VPN dropped (DNS will still resolve). It is
triggered by bursts of connections. `./hpc push` sends everything through a
single tar pipe for exactly this reason; if you do trip it, wait a couple of
minutes and retry.

## Files

| file | role |
| --- | --- |
| `hpc` | the driver — check/setup/push/install/precompile/models/fit/combine/queue/log/ls/pull/cancel/sh |
| `setup-ssh.sh` | per-machine key + `~/.ssh/config` entry |
| `install_pkgs.R` | builds the project R library (`./hpc install`) |
| `precompile.R` | builds the CmdStan precompiled header (`./hpc precompile`) |
| `models.R` | **the model registry** — one entry per model |
| `fit_model.R` | fits the model named by `IGC_MODEL`, one shard per array task |
| `fit.slurm` | array job for the above |
| `combine_model.R` | merges one model's shards, adds `waic`, deletes them |
| `combine.slurm` | job for the above |
| `AGENT.md` | the measurements and the traps — read before changing settings |
| `cogmod_inits_issue.md` | the cold-start init failures and their root cause |
| `hpc.local` | **gitignored** — this machine's account settings, e.g. `IGC_HPC_USER=oc236` |
| `server.md` | **gitignored** — account, keys, OOD URLs |
