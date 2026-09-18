# AGENT.md — running this project on Artemis

Operational notes for driving the Sussex Artemis HPC from this repo. Written
2026-09-17 after getting the fits running end to end; restructured 2026-09-18
into one job per model. `README.md` is the command reference; this file is the
*why*, the traps, and what is still open.

Everything here was measured on the cluster, not assumed.

**Renamed on 2026-09-18**, so that older notes and job logs still make sense:
`make_models.R` -> `fit_model.R` (fits *one* model, named by `IGC_MODEL`),
`make.slurm` -> `fit.slurm`, `combine_models.R` -> `combine_model.R`,
`warmup.R` -> `precompile.R` (it builds the CmdStan precompiled header and has
nothing to do with warm starts or with Stan's warmup), `./hpc submit make` ->
`./hpc fit <model>`. `models.R` is new: the registry every one of them reads.
The warm-start machinery was deleted, per the 2026-09-17 decision to run cold.

---

## 1. Quick start

```
bash analysis/server/setup-ssh.sh   # once per machine, then paste the key into OOD
cd analysis/server
./hpc check                         # VPN + SSH
./hpc install                       # once: build the project R library
./hpc precompile                    # once per toolchain change: build the CmdStan PCH
./hpc push
./hpc fit gam_lnr                   # one array job per model
./hpc queue ; ./hpc log gam_lnr
```

The **GlobalProtect VPN** (portal `bond.sussex.ac.uk`) must be up for any of it.

---

## 2. The working toolchain

These exact pieces work together. Changing any one of them has broken the run
at least once.

| piece | value | why this one |
| --- | --- | --- |
| module | `CmdStanR/0.7.1-foss-2023a-R-4.3.2` | only module stack shipping **mgcv + brms + cmdstanr together**. `R/4.4.1-gfbf-2023b` loads but has neither `mgcv` (needed for `t2()`) nor `cmdstanr`. |
| R | 4.3.2 (GCC 12.3, foss-2023a) | comes with the module |
| brms | 2.21.0 | from the module |
| cmdstanr | **0.9.0**, from the project library | the module's 0.7.1 **cannot read CmdStan 2.39 output** |
| CmdStan | 2.39.0 in `~/.cmdstan` | rebuilt 2026-09-17 for GCC 12.3 |
| project R library | `/mnt/lustre/users/psych/$USER/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3` | holds `datawizard`, `cogmod`, and the newer `cmdstanr`; precedes the module on `R_LIBS` so it shadows it |
| cogmod | **>= 0.3.3** (`dev` branch, commit `e920d44`, installed 2026-09-18) | 0.3.3 fixed the non-finite tail gradient and the init jitter (4.5). `install_pkgs.R` and `fit_model.R` both refuse anything older |

`cogmod` is installed from GitHub and needs an explicit refresh to move:
`./hpc install cogmod`. `IGC_COGMOD_REF` picks the branch — `dev` until 0.3.3
is merged into `main`, at which point change the default in `hpc`.

---

## 3. Traps

Each of these cost a failed run. They are not obvious and they do not announce
themselves clearly.

### 3.1 Compute nodes cannot see the EasyBuild modules

Login nodes have `/mnt/shared/easybuild/modules/all` on `MODULEPATH`; compute
nodes only get `/opt/ohpc/pub/modulefiles`. A job that just does `module load`
dies with *"module(s) exist but cannot be loaded as requested"*, then
`Rscript: command not found`. The `.slurm` files therefore run
`module use "${IGC_EB_MODULES}"` first.

### 3.2 `module` is undefined in SLURM batch scripts

SLURM runs the script non-interactively, where `module` is not a shell
function. The `.slurm` files source `/etc/profile.d/lmod.sh` before any
`module` call.

### 3.3 The CmdStan precompiled header races on a cold array

`~/.cmdstan/<version>/stan/src/stan/model/model_header.hpp.gch/` is a
**directory** holding one ~800 MB variant per compiler/flag combination, e.g.
`model_header_threads_nochecks_12_3.hpp.gch` (GCC 12.3, threading,
`STAN_NO_RANGE_CHECKS`). If the variant we need is absent when an array starts,
every task races to build the same file and all but one die with:

```
error: while reading precompiled header: No such file or directory
```

Worse, the losers leave a **corrupt** variant behind, after which nothing
compiles at all. Recovery is `cmdstanr::rebuild_cmdstan(cores = 8)`.
Prevention is `./hpc precompile`, which builds the variant once on one node.

The variant is keyed on the *flags*, so a `.gch` being present is not the same
as **ours** being present: on 2026-09-18 the directory held the plain and
`nochecks` variants built by the init diagnostics, and none of the threaded one
the fits need. Always look for `model_header_threads_nochecks_12_3.hpp.gch` by
name.

Two non-fixes, both tried:

- `PRECOMPILED_HEADERS=false` in `make/local` — cmdstanr rewrites `make/local`
  from its own `cpp_options` on every compile, discarding it.
- `rm -f .../model/*.gch` — the glob matches the *directory*; `rm -f` refuses
  it ("Is a directory") and the stale variants survive.

Keep `precompile.R`'s `cpp_options` identical to `stan_model_args` in
`fit_model.R`, or a different variant is keyed and the race returns.

### 3.4 Failed tasks used to report COMPLETED

The `.slurm` files ended with an unconditional `echo`, so every task exited 0
regardless of what R did — `sacct` showed `COMPLETED 0:0` for tasks that had
crashed. They now capture `$?` and exit with it. **Do not add a trailing
command after `Rscript` without preserving the status.**

### 3.5 sshd rate-limits connection bursts

Several SSH connections in quick succession get refused with
`kex_exchange_identification: read: Connection reset`. This looks exactly like
a dropped VPN but is not — DNS still resolves. `./hpc push` sends every file
through **one** tar pipe for this reason. If you trip it, wait ~2 minutes.

### 3.6 `file_refit = "never"` reuses stale fits silently

Needed for resumability (see 5.2), but a shard left on disk is reused even if
the formula or data changed. Mitigations, in order:

- `combine_model.R` deletes shards once they are safely combined
- `IGC_FILE_REFIT=always ./hpc fit <model>` forces a clean refit
- delete by hand: `./hpc sh "rm -f <models_dir>/*.rds"`

Probe and test runs **must** use a separate `IGC_MODELS_DIR` so their output can
never be mistaken for production shards. Since 2026-09-18 `./hpc push` also
mirrors rather than merges, so a renamed script cannot be left behind on the
cluster and submitted by accident — the same class of bug one level up.

---

## 4. Measurements

### 4.1 Dataset

| | value |
| --- | --- |
| total rows | 963,885 |
| total participants | 2,221 |
| MullerLyer rows (what is actually fitted) | 323,981 |
| MullerLyer participants | 2,215 |

`fit_model.R` defaults to **all** participants since 2026-09-18 (it used to
default to 30, which risked a production submission silently fitting the test
subset). Override with `IGC_NPARTICIPANTS=30` for a smoke test.

### 4.2 Partition quotas (per user, verified with `sacctmgr`)

| partition | max runtime | max CPUs | max RAM |
| --- | --- | --- | --- |
| `short` | 2 h | 550 | 3.6 TB |
| `general` (default) | 8 h | 400 | 2.7 TB |
| `long` | 8 days | **140** | 900 GB |
| `verylong` | 30 days | 70 | 900 GB |
| `gpu` | 3 days | 300 | 2.1 TB |

Nodes are 128 CPUs / ~478 GB. Defaults if unspecified: 1 CPU per task, 4 GB
per CPU, `general` partition.

**Concurrency is `floor(140 / cpus-per-task)`** — with `--cpus-per-task=16`
that is **8 tasks**, whatever `--array` says, and the cap is per *user*, not per
job: it is shared across everything you have queued. So the array widths of the
jobs running at once should add up to 8. Production is two models x
`--array=1-4`; one model at `--array=1-12` would waste a third of the throughput
running a half-empty second wave.

### 4.3 Memory (1 chain, `save_pars(all = TRUE)` removed)

| participants | MullerLyer rows | MaxRSS |
| --- | --- | --- |
| 30 | 3,841 | 1.69 GB |
| 120 | 15,093 | 1.77 GB |
| 480 | 60,812 | 2.52 GB |

Least-squares fit: **~1.59 GB fixed + ~15 KB/row**, projecting to **~6.5 GB
per chain** at full data. The adapted probes (4.4) give an independent estimate
from different runs — 0.91 GB fixed + 17.6 KB/row, **~6.3 GB per chain** — so
this figure is solid.

Memory is **not** the binding constraint: ~130 chains would fit under the
900 GB cap, well past the 140-CPU ceiling. `--mem=64G` is roughly 5x
over-provisioned for 2 chains, and over-requesting hurts dispatch because TRES
is a priority factor.

Most of the footprint is fixed overhead (R, brms, libraries, the compiled
model), which is why scaling the whole 1.9 GB measured on the toy subset gave a
wildly wrong answer — only the marginal term scales.

Note `sacct` reports `MaxRSS` as the max over *processes*, not the sum — with
several chains per task it shows one chain, not the task total.

### 4.4 Runtime — sublinear; full data is feasible

Measured with **adapted** probes (`warmup=1000` + 100 draws, 2 chains x 8
threads, `gam_lnr` only), which is the production cost regime:

| participants | MullerLyer rows | elapsed | warnings |
| --- | --- | --- | --- |
| 120 | 15,093 | 34:06 | 2% divergences |
| 480 | 60,812 | 1:48:02 | **none** |

Subtracting ~55 s of compile time, 4.03x the rows costs **3.23x** the time —
a scaling exponent of **0.84, i.e. sublinear**. Extrapolated to 323,981 rows
and production's 1300 iterations:

| assumption | per model | all three, sequential |
| --- | --- | --- |
| measured exponent 0.84 | ~8.6 h | **~1.1 days** |
| conservative, linear | ~11.2 h | ~1.4 days |

**Full data fits comfortably inside `long`'s 8-day wall**, under either
assumption, with large margin.

#### The earlier "superlinear" finding was an artifact — do not repeat it

A first pass used `IGC_WARMUP=100` to make probes cheap. Too short for
step-size adaptation, so nearly every transition ran to the max-treedepth
2^10 = 1024 leapfrog steps and the runs measured a far costlier regime:

| run | warmup | max-treedepth hits | divergences | elapsed |
| --- | --- | --- | --- | --- |
| probe 30 | 100 | **100%** | 0 | 2:17 |
| probe 120 | 100 | **98%** | 2% | 6:16 |
| probe 480 (4 chains) | 100 | **75%** | 25% | 2:46:57 |
| adapted 480 (2 chains) | **1000** | **none** | **none** | 1:48:02 |

The unadapted 480 run took **1.5x longer** than the adapted one despite doing
the same 150 vs 1100 iterations — the saturated treedepth dominated everything.
That is what produced a bogus superlinear reading.

**Lesson: shortening warmup is a valid shortcut for measuring memory and a trap
for measuring time.** Always probe runtime at production warmup.

#### 4.4.1 How short can a cold warmup be? (2026-09-17/18, `gam_lnr`, 2 chains x 8 threads, 300 draws)

The trap above is warmup **100**. A separate matrix tested warmup 200-1000,
all cold starts (`cogmod_inits()`, no warm start). One chain per 480-run died
at initialisation (4.5), so those rows are single chains. Times are per chain
from CmdStan's own warmup/sampling split; `s/iter` is the **adapted
sampling-phase** cost; ESS is the median bulk ESS of the population-level
parameters (`b_`, `sd_`, `sds_`), per draw and per 1000 s of total chain time.

| run | participants | warmup | node | warmup s | sampling s | s/iter | leapfrog | ESS/draw | ESS/1000 s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | 120 | 1000 | rtx-02 (zen5, quiet) | 858 | 187 | 0.63 | 255 | 0.31 | 88 |
| S5 | 120 | 500 | a40-16 (shared) | 962 | 349 | 1.17 | 257 | 0.30 | 67 |
| S6 | 120 | 300 | a40-16 (shared) | 669 | 291 | 0.97 | 255 | 0.30 | 95 |
| S8 | 120 | 200 (early windows) | a40-01 (loaded) | 721 | 423 | 1.39 | 255 | 0.28 | 71 |
| S7 | 120 | 300 (early windows) | a40-01 (104/128 busy) | 1414 | 781 | 2.51 | 241 | 0.30 | 41 |
| G1 | 480 | 1000 | a40-03 (106/128 busy) | 9329 | 1876 | 6.25 | 255 | 0.23 | 6.2 |
| G5 | 480 | 500 (2 chains) | a40-10 (most loaded) | 10109 / 12813 | 2553 / 2096 | 7.0-8.5 | 254 | 0.23 | 4.9 |
| G6 | 480 | 300 | a40-16 (shared) | 5059 | 4178 | 13.9 | 511 | 0.35 | 11.4 |
| G8 | 480 | 300 (early windows) | general-02 (128/128) | 2543 | 1363 | 4.54 | 370 | 0.18 | 13.8 |

What this says:

- **Warmup length does not change the quality of a draw.** ESS per draw is
  ~0.30 for warmup 200 to 1000; treedepth, step size and divergence rate
  (0-0.5%) are the same. Warmup 300 is enough for this model at these sizes.
- **Warmup cost is dominated by a fixed phase, not by its length.** The first
  ~100-150 iterations run with the identity metric at max treedepth (1023
  leapfrog steps, ~4x the adapted cost); after the first metric window the
  cost drops to the sampling rate. In sampling-iteration units a cold warmup
  costs roughly `130 x (1023/leapfrog - 1) + warmup`, i.e. ~350 + warmup at
  120 participants. Cutting warmup 1000 -> 300 roughly halves total chain time
  at 480 participants (ESS per second doubles) at no quality cost.
- **Node placement is a 2-3x factor on wall time.** Per-gradient cost was
  2.4 ms on an idle zen5 node and 10 ms on a 104/128-loaded zen3 node for the
  same fit (120 participants). Compare runs by CmdStan's warmup/sampling split
  and leapfrog counts, never by `sacct` elapsed alone.
- **480 participants, adapted: 4.5-14 s per iteration (8 threads)**, 24-27 ms
  per gradient. Leapfrog per iteration varies 255-511 between chains
  (adaptation noise), which puts a 2x band on any single-chain estimate.
- **Within-chain threading scales.** G7 (4 chains x 4 threads) ran at half
  the per-chain speed of the 8-thread runs, so more chains per task buys no
  throughput; 1 chain x 16 threads would halve wall time per chain (but see
  4.5 -- a lone chain that dies at init loses the task).
- **Mixing is limited by a few participants**, not by the population level.
  Max Rhat 1.3-1.9 and min ESS ratio ~0.005 come from one participant's
  `ndt`/`poutlier` random effects trading off (S22 in the 120 subset); the
  population parameters sit at Rhat <= 1.05-1.18 with ESS ~90 per 300 draws.
  Plan for >= 1500 total draws per model to get ESS >= 400 on population terms.

**Full data (2,215 participants, 324k rows) breaks the pattern.** The
single-task pilot with warmup 300 (job 11373507, `tests/F1_cold300_all`,
`general`, 8 h) was read from its CmdStan CSV on the node before it could hit
the wall (copy in `tests/F1_cold300_all/csv/`):

| | value |
| --- | --- |
| warmup 300 iterations | 3 h 22 m (~58-60 s per iteration, 8 threads) |
| adapted step size | **0.0015** (0.010-0.016 at 480) |
| sampling draws | **treedepth 10 / 1023 leapfrog on every draw**, accept 0.94, ~80 s per draw |
| `lp__` over the first 51 draws | rising monotonically 28,499 -> 32,472 |

The rising `lp__` means the chain had **not reached the typical set after 300
warmup iterations**; adaptation ran on transient draws and produced a metric
and step size fit for nowhere. Per gradient it is ~78 ms (1.8x the 480 rows
gave ~1.4x the cost; roughly linear), so the problem is not the gradient but
the length of the cold-start transient, which grows with the number of
participants: ~130 iterations at 120-480, > 300 at 2,215. At ~80 s per
iteration that transient alone is > 7 h. **No warmup short enough to fit an
8 h partition can adapt this model on the full data from a cold start.**

Consequences: full data goes on `long`; warmup must be long enough to get
through the transient *and then* adapt (1000 is the tested value; nothing
shorter has been shown to work at this size); and the 4.4 extrapolation from
480 participants (8.6-11 h per model) is optimistic because the sampler at
full data spends the transient at 1023 leapfrog steps rather than ~255. A
realistic per-chain estimate for warmup 1000 + 300 draws is 7-11 h of
transient plus 4-8 h adapted, i.e. **~12-20 h per model**, with the first
production run as the measurement (watch its `REPORT` line and CmdStan's
warmup/sampling split).

### 4.5 Initialisation is stochastically fragile at scale

At 480 participants with **1 chain** the job died immediately:

```
Chain 1 Rejecting initial value:
Chain 1   Gradient evaluated at the initial value is not finite.
Error: Fitting failed. Unable to retrieve the metadata.
```

The same size with **4 chains** had zero rejections. So it is stochastic, not a
hard break. (Both runs used the short warmup, but initialisation happens before
any adaptation, so this finding is unaffected.) It is fatal rather than a hiccup because brms passes **explicit**
inits from `cogmod_inits()`, and cmdstanr then uses them verbatim — it does not
retry with fresh random draws as it would with default inits.

Measured rate across the 2026-09-17/18 cold-start tests: **7 of 26 chains
(27%) died at initialisation**, at 120, 480 and 2,215 participants alike. A
chain that dies costs nothing but its share of the task; the task continues
with the survivors.

**Fixed in cogmod 0.3.3 (2026-09-18).** The cause was run down to a single
trial per rejection and written up in `cogmod_inits_issue.md`: the `A == 0`
branch of `cogmod_lognormal_acc_ltails()` called Stan's `lognormal_lcdf` /
`lognormal_lccdf`, which are `erfc` alone and underflow ~38 standardized log
units out, giving `-inf` with a `NaN` partial; the outlier mixture kept the
*value* finite, so only the gradient broke. The jittered init only had to put
one trial of ~300k in that region, which is why the rate grew with N. cogmod
0.3.3 routes both through a new `cogmod_log_Phi()` (`erfc` in the body, an
asymptotic tail expansion below `x = -25`) and tiers the init jitter, giving the
`z_*` / `zs_*` / `sd_*` / `sds_*` blocks a fifth of the population jitter and
starting `sds_*` at 0.05 rather than 0.25. The same `NaN` gradient was also
what divergences looked like mid-warmup, so the fix should show up as cleaner
sampling, not only as surviving chains.

**Unverified at production scale** — the first full-data run is the test. Watch
the `.err` files for "Rejecting initial value" and the `REPORT` line's
divergence count.

Practical consequence while it is unverified: **keep 2 chains per task.** One
chain x 16 threads halves the wall per chain and is the obvious next move once
the fix is confirmed at scale, but a single rejection then loses the whole
task.

### 4.6 Sampler health at production settings

At production warmup the sampler is healthy and **gets better, not worse, with
more data**:

| run | warmup | divergences | treedepth |
| --- | --- | --- | --- |
| 30 ppts (smoke, 3 models) | 1000 | 1% (6/600) | none |
| 120 ppts | 1000 | 2% (3/200) | none |
| 480 ppts | 1000 | **0** | none |

The alarming escalation seen earlier (0% -> 2% -> 25% divergences) was entirely
a short-warmup artifact. There is currently **no evidence of a model-geometry
problem**. Re-check at full data, but expect it to be fine.

---

## 5. Design decisions

### 5.1 One job per model, and `models.R` as the single definition

Until 2026-09-18 one array task fitted all three models in sequence and
`IGC_MODELS` filtered them. Now `models.R` is a registry — one entry per model,
holding its illusion and a function returning its `bf()` — and everything else
reads it: `fit_model.R` fits the one named by `IGC_MODEL`, `combine_model.R`
merges the same one, and `./hpc` greps the names out of the file so a typo
fails at submission rather than eight hours in. Adding a model is one entry and
nothing else.

Why one job per model, beyond tidiness: each model gets its own wall clock, a
straggler cannot take the others down with it, the two models run concurrently
under the CPU cap instead of serially inside one task, and per-model job names
and log names (`igc_fit_gam_lnr`, `fit_gam_lnr_<jobid>_<task>.out`) make
`squeue` and the logs readable.

The one contract to keep: `./hpc` reads the model names with
`sed -n 's/^  \([A-Za-z_][A-Za-z0-9_]*\) = list($/\1/p'`, so the declaration
lines in `models.R` must stay at two-space indent in the form `  <name> = list(`.
`./hpc models` prints what it found, which is the check.

### 5.2 `save_pars(all = TRUE)` removed

It retained every Stan parameter including the latents behind roughly
6 x n_participants group-level coefficients. Only moment-matched `loo` needs
it; this project uses `waic`. Restore it if `loo(moment_match = TRUE)` is ever
wanted.

### 5.3 `file_refit = "never"`

Stan cannot checkpoint mid-chain, but the shard is a good checkpoint: each
array task writes its own `.rds` as it completes. With `"never"`, a resubmitted
array skips shards that finished and resumes. With `"always"` (the old value) a
job killed at the wall restarted from nothing. See 3.6 for the staleness
trade-off; brms also offers `"on_change"`, which refits automatically when
formula, data or prior change.

### 5.4 The LNR formula must fix `sigmabias`

`cogmod_lnr()` gained a `sigmabias` dpar. Omitting it from the `bf()` does
**not** error — it silently becomes a freely estimated start-point range with a
`lognormal(-0.35, 0.75)` prior. cogmod's own comment: *"0 is the plain LNR."*
The formula therefore sets `sigmabias = 0` explicitly. Watch for this class of
bug whenever cogmod adds a parameter — and note that as of 0.3.3 the
`sigmabias = 0` path is also the cheaper one (the tails are computed singly
rather than as a discarded pair, ~20% less per gradient).

### 5.5 No warm start; the machinery is gone

Decision 2026-09-17, implemented 2026-09-18: production fits are cold starts
from `cogmod_inits()` with Stan's own adaptation. `fit_model.R` no longer reads
`warmstart.csv`, and `IGC_WARMSTART` / `IGC_INIT_BUFFER` / `IGC_WINDOW` /
`IGC_TERM_BUFFER` / `IGC_ADAPT_ENGAGED` are gone with it — they existed to make
a supplied inverse metric survive adaptation, which is meaningless without a
warm start. `cogmod_warmstart()` still exists in cogmod if the decision is ever
revisited; `analysis/warmstart.csv` and the chunk in
`analysis/1_modelcomparison.qmd` that writes it are untouched.

### 5.6 Everything is parameterised by environment variable

No user paths are hard-coded; `./hpc` passes them via `sbatch --export`.

| variable | default | purpose |
| --- | --- | --- |
| `IGC_HPC_USER` | `dmm56` | cluster account |
| `IGC_USERS_DIR` | `/mnt/lustre/users/psych/$USER/IGComputational` | code + models |
| `IGC_SCRATCH_DIR` | `/mnt/lustre/scratch/psych/$USER/IGComputational` | job logs |
| `IGC_MODELS_DIR` | `$IGC_USERS_DIR/models` | fitted shards — **give test runs their own** |
| `IGC_MODEL` | none | which model; set by `./hpc fit` / `./hpc combine` |
| `IGC_R_MODULE` | `CmdStanR/0.7.1-foss-2023a-R-4.3.2` | module to load |
| `IGC_R_LIBS` | project library (4.3) | extra packages |
| `IGC_EB_MODULES` | `/mnt/shared/easybuild/modules/all` | EasyBuild tree |
| `IGC_COGMOD_REF` | `dev` | branch/tag `./hpc install` tracks |
| `IGC_FILE_REFIT` | `never` | `always` forces a clean refit |
| `IGC_KEEP_SHARDS` | unset | set to keep shards after combining |
| `IGC_NPARTICIPANTS` | `all` | subset size for tests, e.g. `30` |
| `IGC_WARMUP` | `1000` | warmup iterations |
| `IGC_SAMPLES` | `500` | post-warmup draws per chain |
| `IGC_CHAINS` | `2` | chains per array task; threads per chain is `cpus / chains` |

Every fit prints one `REPORT <model> | wall | n_leapfrog | treedepth | stepsize |
divergent | accept | max Rhat | min neff_ratio` line to the `.out` log, so runs
can be compared without pulling the `.rds`.

---

## 6. The production configuration

Settled 2026-09-18 from the measurements above. `fit.slurm`'s `#SBATCH`
defaults *are* this configuration, so the production run is two bare commands:

```bash
./hpc push
./hpc fit gam_lnr
./hpc fit gam_lnr6
./hpc fit gam_ddm4
```

| setting | value | why |
| --- | --- | --- |
| partition | `long`, `--time=2-00:00:00` | 12-20 h per chain expected (4.4.1); 8 h is out of reach and there is no 24 h tier |
| participants | `IGC_NPARTICIPANTS=all` (default) | 2,215 participants, 323,981 rows |
| warmup | `IGC_WARMUP=1000` | the only value shown to adapt at full data; 300 provably does not (4.4.1) |
| draws | `IGC_SAMPLES=500` | 8 chains x 500 = 4,000 draws per model; at ESS/draw ~0.3 that is ~1,200 ESS on population terms |
| array | `--array=1-4` per model | 4 x 16 = 64 CPUs per model against `long`'s 140-CPU per-user cap, so 8 of the 12 production tasks run at once and the rest start as they finish (`--time` is per task, so pending costs nothing) |
| chains | `IGC_CHAINS=2`, 8 threads each | threading scales nearly linearly, so more chains per task buys nothing; one chain per task risks the whole task on one init (4.5) |
| memory | `--mem=32G` | ~6.5 GB per chain at full data (4.3), so ~13 GB used; the rest is headroom against an OOM killing a 20 h job |

`gam_ddm5` is defined in `models.R` and submittable, but is not part of the
final run. It carries one more 2-D smooth than DDM-4, so give it its own
`--time` if it is ever wanted.

`gam_lnr6` (added 2026-09-18) is the LNR with `sigmabias` estimated and
smoothed like the other parameters. It is the one model expected to be *slower*
per gradient than its sibling: cogmod 0.3.3 made the `sigmabias = 0` path ~20%
cheaper by computing a single tail, and the free-`sigmabias` path ~15% dearer
in exchange for a corrected gradient. If anything straggles, expect it to be
this one.

### 6.1 What to watch on the first production run

It doubles as the measurement that closes open question 7.1. In the `.out`:

- `n_leapfrog` on the `REPORT` line **well below 1023**, and a step size around
  0.01. Still at 1023 means warmup 1000 did not clear the cold-start transient
  either, and the answer is better inits or a higher `max_treedepth`, not more
  hours.
- `divergent` near 0. At production warmup the sampler was clean at every size
  tested (4.6), and cogmod 0.3.3 removed the `NaN` gradient that made some of
  the earlier divergences (4.5).
- In the `.err`: no "Rejecting initial value". Any at all means the 0.3.3 fix
  is incomplete at this scale — reopen `cogmod_inits_issue.md`.
- `max Rhat` and `min neff_ratio` are expected to look bad (1.3-1.9, ~0.005)
  because one or two participants' `ndt`/`poutlier` trade off; check the
  population-level terms rather than the global maximum (4.4.1).

### 6.2 Levers if it is too slow

- **More draws are the cheapest samples available** once warmup is paid. With
  warmup 1000 (77% overhead) 300 -> 1000 draws roughly triples the yield for
  ~1.5x the wall. On `long` lengthen the chain rather than adding tasks.
- **1 chain x 16 threads halves the wall per chain**, and becomes reasonable
  once the 0.3.3 init fix is confirmed at scale (4.5, 6.1).
- **Not** a shorter warmup: at full data the cold-start transient alone exceeds
  300 iterations at ~80 s each (4.4.1).
- **Not** `verylong`: 30 days at half the CPUs is a net loss.

---

## 7. Open questions

1. **How long is the full-data cold-start transient, and does warmup 1000
   clear it?** The warmup-300 pilot (4.4.1) shows > 300 iterations at 1023
   leapfrog steps (~80 s each) and no adaptation. Warmup 1000 has never been
   run at full data. The first production `gam_lnr` job answers this — see 6.1.
2. **Does cogmod 0.3.3 hold at production scale?** (4.5) The init rejections
   and their root cause are understood and fixed, and a 30-participant smoke
   test on 2026-09-18 started all four chains cleanly, but the failure rate
   grew with N, so the full-data run is the real test. If it holds, 1 chain x
   16 threads per task becomes the better layout.
3. **Do `gam_ddm4` / `gam_ddm5` cost the same per iteration as `gam_lnr`?** At
   30 participants they did (DDM-5 ~1.4x). Not yet measured at scale — the two
   production jobs give the DDM-4 half of the answer.
4. **Does `gam_lnr6` pay for itself?** A free `sigmabias` costs ~15% more per
   gradient and adds a sixth 2-D smooth; whether the start-point range is
   identified at all on this data is what the `waic` comparison against
   `gam_lnr` is for.

*Answered:* `combine_model.R` — run against real shards on 2026-09-18 and
now working; it turned up three bugs the old `combine_models.R` had never had
the chance to show (no `library(cogmod)`, so `add_criterion("waic")` could not
resolve `log_lik_cogmod_lnr()`; a hardcoded `ndraws = 1500` that errors on any
fit with fewer draws; and the inherited `$file` slot that made
`add_criterion()` overwrite shard 1 with the combined fit). Runtime scaling
(4.4, sublinear), memory (4.3, ~6.3-6.5 GB/chain),
initialisation (4.5, root cause found and fixed in cogmod 0.3.3), sampler
health (4.6, clean), warmup length (4.4.1: 300 is enough below full data, and
not enough at it). *Closed by decision:* warm start — not used (2026-09-17),
machinery removed (2026-09-18); `verylong` — not worth it.

## 8. File map

| file | role |
| --- | --- |
| `hpc` | the driver — check/setup/push/install/precompile/models/fit/combine/queue/log/ls/pull/cancel/sh |
| `setup-ssh.sh` | per-machine key + `~/.ssh/config` entry |
| `install_pkgs.R` | builds the project R library (`./hpc install`); enforces cogmod >= 0.3.3 |
| `precompile.R` | builds the CmdStan precompiled header (`./hpc precompile`) |
| `models.R` | **the model registry** — one entry per model, read by everything else |
| `fit_model.R` | fits the model named by `IGC_MODEL`, one shard per array task |
| `fit.slurm` | array job for the above; its `#SBATCH` defaults are the production config |
| `combine_model.R` | merges one model's shards, adds `waic`, deletes them |
| `combine.slurm` | job for the above |
| `README.md` | command reference, partition quotas, PCH detail |
| `cogmod_inits_issue.md` | the cold-start init failures, root cause, and the 0.3.3 fix |
| `server.md` | **gitignored** — account, keys, OOD URLs |
