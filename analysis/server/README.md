# Running the models on Artemis (Sussex HPC)

> **Cluster-level instructions live in the lab HPC hub**:
> <https://github.com/RealityBending/Lab/tree/main/hpc> (start at `hpc/README.md`;
> on Dom's machines: `~/Dropbox/RealityBendingLab/Lab/hpc/`).
> Prefer a local clone of `RealityBending/Lab` if there is one — it also holds
> your gitignored `hpc/private/` notes. The hub covers access, storage,
> partitions and quotas, the R/Stan toolchain, job conventions, troubleshooting
> and housekeeping, and **wins over anything here that contradicts it**; this
> file should only hold what is specific to this project.

Model fitting is too heavy for a laptop, so it runs as SLURM array jobs on
**Artemis**. The `./hpc` script in this folder wraps the whole loop — push
code, submit, watch, pull results — over SSH.

**One job fits one model.** `models.R` is the registry of what a model *is*;
`./hpc fit <model>` submits an array for it, each array task writing one shard;
`./hpc combine <model>` merges that model's shards into a single fit.

## The cluster, and one-time setup

Access (VPN, Open OnDemand, SSH), storage, partitions and quotas:
[hub `artemis.md`](https://github.com/RealityBending/Lab/blob/main/hpc/artemis.md). Setting up a machine (SSH key,
`artemis` alias) or a new account (R library, CmdStan, precompiled header):
[hub `setup.md`](https://github.com/RealityBending/Lab/blob/main/hpc/setup.md).

This project's directories on the cluster:

| path | use |
| --- | --- |
| `/mnt/lustre/users/<group>/<user>/IGComputational/` | code, `models/` (shards), `models/combined/` |
| `/mnt/lustre/scratch/<group>/<user>/IGComputational/` | job logs (`fit_<model>_<job>_<task>.out`) |

Once per account, for this project:

```bash
cd analysis/server
./hpc check && ./hpc setup && ./hpc install && ./hpc precompile
```

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
./hpc combine gam_lnr      # merge its shards and add loo (shards are kept)
./hpc pull                 # bring combined fits into analysis/models/
```

`push` **mirrors** rather than merges: a top-level `*.R` / `*.slurm` on the
cluster that no longer exists locally is deleted, so a renamed script can never
be submitted by accident. `models/`, `tests/` and the logs are untouched.

### How to check job status, cheaply

Report status the hub's way ([`jobs.md#status-reports`](https://github.com/RealityBending/Lab/blob/main/hpc/jobs.md#status-reports)):
running since when, a per-shard table of each chain's iterations, an ETA
against that model's warmup + samples, and a warning if a shard risks the
wall. This driver has no `progress` command yet, so:

```bash
./hpc queue                # one line per array task: STATE, TIME, TIME_LEFT, reason if PENDING
./hpc sh "grep -h 'Iteration\|REPORT' \$IGC_SCRATCH_DIR/fit_<model>_*.out | tail -40"
```

The `REPORT ...` line is printed once at the end of a shard (see `AGENT.md`
§6.1 for what to expect in it). Seeing only `Iteration: 1 / 1500` for hours is
normal here: at full data warmup 1000 alone costs 12-17 h (`AGENT.md` §4.4.1,
§7 Q1). Trust `TIME` in `./hpc queue` over the absence of a fresh line.

Extra arguments to `fit` and `combine` are passed to `sbatch`, so a smoke test
is:

```bash
IGC_MODELS_DIR=/mnt/lustre/users/psych/dmm56/IGComputational/smoke \
IGC_NPARTICIPANTS=30 IGC_WARMUP=300 IGC_SAMPLES=100 \
  ./hpc fit gam_lnr --array=1-2 --partition=short --cpus-per-task=8 --mem=16G
```

`--partition=short` is the only wall-clock setting needed: with no `--time`,
the task gets that partition's maximum (2 h), which is all a smoke test can
use anyway.

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
`--cpus-per-task=16`, `--mem=32G`, `--partition=long` and no `--time` (so each
task gets `long`'s 8-day maximum), with `IGC_NPARTICIPANTS=all`,
`IGC_WARMUP=1000`, `IGC_SAMPLES=500` and `IGC_CHAINS=2`. That is 4 shards x 2
chains x 500 = **4,000 draws per model**.

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
CPU: `--mem=128G` on `general`, because `add_criterion()` builds a
draws x 323,981 pointwise log-likelihood matrix (7.2 GB over all draws) through
an R-level `log_lik` called once per response. Measured peak was 45 GB.

Since 2026-09-20 it adds **`loo` over every draw** and **keeps the shards**.
Both were measured rather than assumed: `loo` over all 3,000 draws of the
3-shard `gam_lnr` took 9.5 minutes against an 8 h wall, and costs only 1.9x
what 500 draws cost, so subsampling trades real precision for almost no time.
`cores` does not help (6% at 16 CPUs). AGENT.md §4.8 has the table.

Override per run when a model actually overruns — not pre-emptively:

```bash
IGC_CRITERION=waic ./hpc combine gam_ddm7          # cheaper criterion
IGC_CRITERION_NDRAWS=1500 ./hpc combine gam_ddm7   # or fewer draws
IGC_CRITERION=none ./hpc combine gam_ddm7          # merge only
```

Keeping shards reopens the trap that deleting them used to close: `file_refit =
"never"` means a shard on disk is silently reused even when the formula or the
data changed. **A changed parametrisation needs its own `IGC_MODELS_DIR`, or
`IGC_FILE_REFIT=always`.**

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

| model | family | distributional parameters with a 2-D smooth | fixed | |
| --- | --- | --- | --- | --- |
| `gam_lnr` | `cogmod_lnr()` | drift, `nuone`, `sigmazero`, `sigmaone`, `ndt` | `sigmabias = 0` |
| `gam_lnr6` | `cogmod_lnr()` | as LNR plus `sigmabias` | — |
| `gam_ddm4` | `cogmod_ddm()` | drift, `boundary`, `bias`, `ndt` | `sigmadrift`/`sigmabias`/`sigmandt` = 0 |
| `gam_ddm5` | `cogmod_ddm()` | as DDM-4 plus `sigmadrift` | `sigmabias`/`sigmandt` = 0 |
| `gam_ddm7` | `cogmod_ddm()` | all seven: as DDM-5 plus `sigmabias`, `sigmandt` | — | ⚠ **not viable at full data — subsample only** |
| `gam_rdm` | `cogmod_rdm()` | drift, `driftone`, `boundary`, `ndt` | `sigmabias = 0` |
| `gam_rdm5` | `cogmod_rdm()` | all five: as RDM plus `sigmabias` | — |
| `gam_lba` | `cogmod_lba2()` | drift, `driftone`, `sigmaone`, `sigmabias`, `boundary`, `ndt` | `sigmazero = 1` |

Every smooth is `t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs =
c("cr", "cr")) + (1 | Participant)`; `poutlier` is `1 + (1 | Participant)`.
Two modelling choices in that table are worth knowing about before reading any
output:

- **`gam_lba` is `cogmod_lba2`, not `lba1`.** `lba1` has a single drift and no
  second accumulator, so it cannot model the choice that `dec(Error)` carries.
- **`sigmazero = 1` in `gam_lba` is a scaling constraint, not a simplification.**
  The LBA is identified only up to a scale — multiply every drift, the
  start-point range and the boundary by one constant and the likelihood does
  not move — so one parameter must be pinned, and cogmod's own convention is
  the first accumulator's drift SD. `sigmaone` stays free, which is what lets
  the accumulators differ. Freeing `sigmazero` without substituting another
  constraint leaves a ridge for the chains to wander along.

`gam_rdm` fixes `sigmabias = 0`, the plain racing diffusion with both
accumulators starting from the same point; `gam_rdm5` frees it, and stands to
`gam_rdm` as `gam_lnr6` does to `gam_lnr`. Unlike the LBA the RDM needs no
scaling constraint — its diffusion coefficient is fixed internally — so freeing
`sigmabias` there does not open the ridge that freeing `sigmazero` would in
`gam_lba`.

### `gam_ddm7`: do not submit it at full data

Smoke-tested on 2026-09-18 and **not viable at full data**. At 30 participants
it did not reach iteration 100 in 50 minutes, where every other model in the
registry finished all 400 in 6-11 minutes — upwards of 25x the per-iteration
cost on 1/74th of the production data, which scales to roughly a fortnight per
chain at full data, past `long`'s 8-day ceiling.

The cause was measured later the same day and it is **not** the geometry
problem this section first guessed at. `cogmod_ddm_decision_lpdf()` routes on an
exact-zero test, and estimating `sigmabias` or `sigmandt` — as opposed to fixing
either at `0` in `bf()` — leaves the analytic Wiener density for Stan's adaptive
numerical quadrature, at 18x (one freed) or 55x (both) the cost of the path
`gam_ddm5` takes. `sigmadrift` is analytic and costs 2.8x, which is why
`gam_ddm5` is cheap.

What follows from that, before anyone spends a run rediscovering it:

- **Tighter priors do not help.** Driving both to `1e-5` still costs 26x, and
  that is a floor — the fast path tests for *exact* zero.
- **The smooths are not the problem.** The cost is per observation, so making
  the variabilities intercept-only changes nothing.
- **A warm start does not help either**, at full data: it buys back warmup,
  and the retained draws alone are ~18 days.
- **Subsampling does.** `gam_ddm7` is ~2.7-4.7 days per chain at 200
  participants with 1 chain x 16 threads, which is the way to get a look at it.
  Give it its own `IGC_MODELS_DIR`.

`gam_ddm4` and `gam_ddm5` are unaffected. Full numbers and the proposed cogmod
fix in `AGENT.md` §4.7.1 and `cogmod_ddm_cost_issue.md`.

Adding a model is one entry in `models.R` and nothing else. Which account is
fitting what is deliberately not recorded here — it changes, and `./hpc queue`
answers it live.

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

A second account is no longer the *only* way to get a second allowance:
`dmm56` can also use the `sussexneuro` partition, whose quota is independent of
`long`'s (`AGENT.md` §4.2.1; rules in the hub's `artemis.md#sussexneuro`). The
two stack, so the fastest arrangement is a colleague on `long` plus a
`--partition=sussexneuro` job here.

Their setup, once:

```bash
bash <Lab>/hpc/scripts/setup-ssh.sh oc236             # their key, their ssh alias (lab hub)
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
# one account
./hpc fit gam_lnr ; ./hpc fit gam_lnr6 ; ./hpc fit gam_ddm4
# the other
./hpc fit gam_rdm ; ./hpc fit gam_rdm5 ; ./hpc fit gam_lba
```

Three jobs x 4 tasks x 16 CPUs is 192 against the 140-CPU cap, so the third
waits for the first to finish. That is fine and costs
nothing (`--time` is per task), but if the wall matters, submit two and hold the
third, drop to `--array=1-2` and take 2,000 draws per model instead of 4,000,
or send the third to `--partition=sussexneuro`, which does not draw on that cap
at all.

Both of those three smoke-tested clean at 30 participants (`AGENT.md` §4.7):
`gam_rdm` is the best-behaved model in the registry, and `gam_lba`'s 6%
divergences at warmup 300 are the one number to check in its production
`REPORT`, where warmup 1000 is expected to clear them.

A model that only one of you fits still has to be **defined in the shared
`models.R` and committed**, or `./hpc fit` will reject the name. That is the
point: the definition is reviewed and version-controlled, the account is not.

### Collecting the results

No account can read another's directories, so the files are handed over
([hub `jobs.md#sharing-results-between-accounts`](https://github.com/RealityBending/Lab/blob/main/hpc/jobs.md#sharing-results-between-accounts)).
Whoever fitted the model runs:

```bash
./hpc combine <model>   # adds loo; the shards stay put
./hpc pull              # combined/*.rds -> their own analysis/models/
```

and then sends the file. Pull `combined/*.rds`, not raw shards — `combine`
keeps the shards, so `models/` holds both and each shard is ~214 MB. A
combined fit is large: the 3-shard `gam_lnr` was 428 MB, so a 4-shard one with
`loo` attached is ~600 MB.

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
| `IGC_CRITERION` | `loo` | `waic` is cheaper; `none` skips it |
| `IGC_CRITERION_NDRAWS` | all | subsample the draws the criterion uses |
| `IGC_DELETE_SHARDS` | unset | set to `1` to drop shards after combining |
| `IGC_COGMOD_REF` | `dev` | which cogmod branch/tag `./hpc install` tracks |

Each fit prints a `REPORT ...` line (wall time, mean leapfrog steps,
treedepth, step size, divergences, max Rhat, min ESS ratio) to its `.out`, so
runs can be compared without pulling the `.rds`.

## Precompile the CmdStan header before a cold array

Run `./hpc precompile` once after any change to the toolchain, the CmdStan
version, or the `stan_model_args` in `fit_model.R`, and keep `precompile.R`'s
`cpp_options` identical to them. Our fits need the
`model_header_threads_nochecks_12_3.hpp.gch` variant (GCC 12.3,
`threading()`, `STAN_NO_RANGE_CHECKS`). Why, and how to recover from the race
(it killed 10 of 12 tasks in job 11366858 on 2026-09-17):
[hub `toolchain.md#precompiled-header`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md#precompiled-header).

## Partitions and resource limits

Quotas and partition rules are in [hub `artemis.md`](https://github.com/RealityBending/Lab/blob/main/hpc/artemis.md). For
this project:

- **Concurrency on `long` is `floor(140 / 16) = 8` tasks** across *all* the
  account's jobs there — which is why production is 4 tasks per model for two
  models rather than 8 for one. Anything beyond sits in
  `PENDING (QOSMaxCpuPerUserLimit)`.
- Memory is not the binding constraint: 8 x 32 GB = 256 GB against `long`'s
  900 GB.
- **`fit.slurm` sets no `--time`**, so each task gets `long`'s 8 days; it
  shipped with `--time=2-00:00:00` and real shards then measured 42.7-44 h
  (`AGENT.md` §3.7).
- **`sussexneuro`** runs a further model alongside the `long` arrays on a
  separate quota (`AGENT.md` §4.2.1), e.g. `./hpc fit gam_rdm
  --partition=sussexneuro`. Its 256 CPUs are a group pool — look before
  sizing (hub `artemis.md#sussexneuro`).
- Check where you stand: `./hpc sh "squeue --Format=JobID,State,Reason,PriorityLong,Partition -u dmm56"`.

## R environment on the cluster

The jobs load **one module**, `CmdStanR/0.7.1-foss-2023a-R-4.3.2` — the only
Artemis stack shipping `mgcv` (for `t2()`), `brms` and `cmdstanr` together —
plus the account's project library for what it lacks (`datawizard`, `cogmod`,
a newer `cmdstanr`), built by `install_pkgs.R` via `./hpc install`:

```
/mnt/lustre/users/psych/$IGC_HPC_USER/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3
```

The library is shared with every other project on the account (FakeArt
included). Stack, versions, module quirks and the brms trap on fresh accounts:
[hub `toolchain.md`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md).

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

Data is read straight from the GitHub raw URLs in `fit_model.R` — the compute
nodes do have outbound internet, so nothing needs pushing but code.

`.gitattributes` forces LF on `*.slurm` so a CRLF file never reaches the
cluster (which fails with `bad interpreter`). `./hpc push` strips CR as well,
belt and braces.

## Files

| file | role |
| --- | --- |
| `hpc` | the driver — check/setup/push/install/precompile/models/fit/combine/queue/log/ls/pull/cancel/sh |
| `install_pkgs.R` | builds the project R library (`./hpc install`) |
| `precompile.R` | builds the CmdStan precompiled header (`./hpc precompile`) |
| `models.R` | **the model registry** — one entry per model |
| `fit_model.R` | fits the model named by `IGC_MODEL`, one shard per array task |
| `fit.slurm` | array job for the above |
| `combine_model.R` | merges one model's shards, adds `loo`, keeps them |
| `combine.slurm` | job for the above |
| `AGENT.md` | the measurements and the traps — read before changing settings |
| `cogmod_inits_issue.md` | the cold-start init failures and their root cause |
| `cogmod_ddm_cost_issue.md` | why `gam_ddm7` is 55x dearer per gradient, and the cogmod fix for it |
| `hpc.local` | **gitignored** — this machine's account settings, e.g. `IGC_HPC_USER=oc236` |
| `server.md` | **gitignored** — local notes on this project's cluster dirs |
