# AGENT.md — running this project on Artemis

> **Cluster-level instructions live in the lab HPC hub**:
> <https://github.com/RealityBending/Lab/tree/main/hpc> (start at `hpc/README.md`;
> on Dom's machines: `~/Dropbox/RealityBendingLab/Lab/hpc/`).
> Prefer a local clone of `RealityBending/Lab` if there is one — it also holds
> your gitignored `hpc/private/` notes. The hub covers access, storage,
> partitions and quotas, the R/Stan toolchain, job conventions, troubleshooting
> and housekeeping, and **wins over anything here that contradicts it**; this
> file should only hold what is specific to this project.

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

**Renamed on 2026-10-13: the MullerLyer models carry `_mullerlyer`** like the
other illusions (`gam_lnr` -> `gam_lnr_mullerlyer`, fit file
`gam_lnr_mullerlyer_MullerLyer.rds`). The formulas moved to `igc_formulas` in
`models.R`, keyed by the old bare names, and are unchanged. Everything in this
file written before that date uses the old names: read `gam_lnr` there as
"the LNR on MullerLyer" (or the LNR family, as the context says).

---

## 1. Quick start

```
# once per machine: the lab hub's setup.md (scripts/setup-ssh.sh)
cd analysis/server
./hpc check                         # VPN + SSH
./hpc install                       # once: build the project R library
./hpc precompile                    # once per toolchain change: build the CmdStan PCH
./hpc push
./hpc fit gam_lnr_mullerlyer                   # one array job per model
./hpc queue ; ./hpc log gam_lnr_mullerlyer
```

---

## 2. The toolchain

The cluster stack — module, R, brms, cmdstanr, CmdStan, the project R library
shared by every project on the account — and why each piece is what it is:
[hub `toolchain.md`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md). What is specific to this project:

- The module must be `CmdStanR/0.7.1-foss-2023a-R-4.3.2`: the `t2()` smooths
  need **mgcv**, which `R/4.4.1-gfbf-2023b` does not ship.
- **cogmod 0.3.3** (installed 2026-09-18 and refreshed to `d04c7f8` on
  2026-09-20) fixed the non-finite tail gradient and the init jitter (4.5).
  Every production fit so far was made with it or with something older. It
  was merged into `main` on 2026-09-26.
- **New fits use cogmod 0.3.4** (`dev`, `e8c0555`, 2026-10-02). That is the
  user's decision of 2026-10-02. Once every wanted model is fitted, the user
  will decide whether to refit the ones made with older versions. 0.3.4
  changes the priors of every model with a smooth. In cogmod's `NEWS.md` it:
  - gives `sds` an `exponential(rate)` prior, with the rate set per smooth from
    its own basis, so that the prior is `exponential(1)` in link units;
  - gives the unpenalised `bs_*` columns on a family-named dpar that dpar's
    slope prior, in link units;
  - in the race families, gives `mu`'s participant SD the same
    `exponential(1)` as the other accumulator's;
  - starts smooth SDs at 0.1 link units of wiggle and LBA2's `driftone` at 1;
  - scales the init jitter on `bs_*` blocks by how far their design columns
    reach.

  How that is enforced, from 2026-10-02:
  - `IGC_COGMOD_REF` stays `dev`, so `./hpc install cogmod` installs whatever
    `dev` holds: 0.3.4, then 0.3.5 from 2026-10-12 (next bullet).
  - **The floor is 0.3.5** (0.3.4 from 2026-10-02 to 2026-10-12).
    `install_pkgs.R` and `fit_model.R` refuse anything older.
    `IGC_COGMOD_MIN=0.3.4 ./hpc fit <model>` lets a shard of an older model
    run while the library is still 0.3.4; it cannot bring 0.3.4 back after
    0.3.5 is installed.
  - **Each new shard records its cogmod** as `m$cogmod` (version and commit),
    which `combine_model.R` carries into the combined fit. A fit without
    `$cogmod` predates the stamp, so it is 0.3.3 or older. That is how to list
    the refit candidates later. Fits submitted before 2026-09-20 also predate
    0.3.3's `sds` fix (5.7).
  - **`combine_model.R` refuses shards of mixed cogmod versions**, an
    unstamped shard counting as one version. `brms::combine_models()` does not
    compare priors, so without this check a mixed set would combine silently.

  The library is shared by all of an account's jobs, FakeArt's included.
  Installing a new cogmod while an array still has tasks pending gives those
  shards other priors than their siblings. The combine check now catches that, but it
  costs the shard. The order is:
  1. wait until no array on the account has pending tasks;
  2. `./hpc push`;
  3. `./hpc install cogmod`;
  4. then submit.

  After a push, a pending task of an older array meets the floor and stops in
  its first second. That is the intended failure.

  A `loo` comparison between a 0.3.4 fit and an older one compares models on
  different priors. At 324k rows the likelihood dominates, but smooth `sds`
  and the `bs_*` slopes moved a lot (5.7).
- **New fits use cogmod 0.3.5** (`dev`, user's decision of 2026-10-12),
  which implements part of `cogmod_priors_review.md`. cogmod refitted each
  proposal on its vignettes' data before implementing it
  (`benchmarks/driftone_priors/README.md` in cogmod). In cogmod's `NEWS.md`
  it:
  - centres `poutlier` at 0.25% rather than 0.7%: `normal(-6, 1)` on the
    logit link, `exponential(280)` when `poutlier` is left out of `bf()`.
    Every model's priors change;
  - gives the DDM's `boundary` `normal(1, 1)` on its intercept and
    `normal(0, 0.5)` on its slopes, and the `bias` slopes `normal(0, 0.5)`.
    Both participant SDs get `exponential(1)` instead of brms's
    `student_t(3, 0, 2.5)`, so `fit_model.R`'s blanket `b` prior now reaches
    only `mu`;
  - does **not** re-centre the RDM's and LBA's `driftone` priors (R2-R4 of the
    review). Refitted, those predicted the same error rates and RTs, cost
    9-26% more leapfrog steps per draw on the RDM, and added divergences on
    the LBA2. A `driftone` posterior in its prior's lower tail is expected for
    the error accumulator (cogmod's `?cogmod_priors`), not a flag to act on.

  The order above applies, with one step in front: cogmod's 0.3.5 has to be
  committed and pushed to `dev` before `./hpc install cogmod` can install it.

  **Installed on `dmm56` on 2026-10-13: 0.3.5 at `eb4c6a3`** ("LNR inits",
  the `dev` head that day). Besides the priors above, it carries the LBAs'
  density and gradient fixes and starts the LNR's error accumulator slower
  (`nuone = 0` against `mu = 0.7`). Nothing was pending. The one job running
  on the account, a cogmod `stanli` benchmark that read cogmod from the same
  library, was cancelled first at the user's request. Until then `./hpc
  install` passed a floor of 0.3.4 to `install_pkgs.R` whatever the file
  said, so its final check could not catch a stale 0.3.4. It passes 0.3.5
  now.
- Keep `precompile.R`'s `cpp_options` identical to `stan_model_args` in
  `fit_model.R`, or a different precompiled-header variant is keyed and the
  race comes back ([hub](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md#precompiled-header)).

---

## 3. Traps

Only 3.6 and 3.9 are specific to this project. The rest are cluster facts and moved to
the lab hub on 2026-09-22; the numbers are kept so older references resolve.

| § | trap | now in the hub |
| --- | --- | --- |
| 3.1, 3.2 | compute nodes cannot see the EasyBuild modules; `module` undefined in batch scripts | [`toolchain.md#modules`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md#modules) |
| 3.3 | the CmdStan precompiled header races on a cold array | [`toolchain.md#precompiled-header`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md#precompiled-header) |
| 3.4 | failed tasks reported COMPLETED | [`jobs.md#the-slurm-script`](https://github.com/RealityBending/Lab/blob/main/hpc/jobs.md#the-slurm-script) |
| 3.5 | sshd rate-limits connection bursts | [`troubleshooting.md#ssh-connection-reset`](https://github.com/RealityBending/Lab/blob/main/hpc/troubleshooting.md#ssh-connection-reset) |
| 3.7 | do not set `--time` below the partition's max for a job of uncertain length (this project's 42.7-44 h `gam_lnr` shards against a 2-day `--time`, 2026-09-20, are the example there) | [`jobs.md#wall-time`](https://github.com/RealityBending/Lab/blob/main/hpc/jobs.md#wall-time) |
| 3.8 | a second account's fits fail to compile: `install_pkgs.R` does not reproduce the working brms, which `dmm56` borrows through `R_LIBS_USER` (found 2026-09-21; pinning brms in `install_pkgs.R` still to do) | [`toolchain.md#the-brms-trap`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md#the-brms-trap) |

### 3.6 `file_refit = "never"` reuses stale fits silently

Needed for resumability (see 5.2), but a shard left on disk is reused even if
the formula or data changed. Mitigations, in order:

- **`fit_model.R` reuses a shard only if it belongs to this run** (since
  2026-10-13). It must have a cogmod stamp at or above the floor, and the
  run's warmup and draw count. Anything else stops the task in its first
  minutes ("is from another run"). It does not refit over the shard, which
  may be the only copy. A different formula or dataset under the same cogmod
  and run length still gets through, so the next point stands.
- **`./hpc archive`** before a new generation of fits moves every shard and
  combined, prediction and priors file into `models/old_<date>/`, out of
  reach of every script (README → "Starting a new generation of fits").
  `dmm56`'s pre-refit generation is in `models/old_20261013/`.
- `IGC_FILE_REFIT=always ./hpc fit <model>` forces a clean refit
- `combine_model.R` no longer deletes shards after combining (since
  2026-09-20, 4.8). `IGC_DELETE_SHARDS=1` brings that back.

Probe and test runs **must** use a separate `IGC_MODELS_DIR` so their output can
never be mistaken for production shards. Since 2026-09-18 `./hpc push` also
mirrors rather than merges, so a renamed script cannot be left behind on the
cluster and submitted by accident — the same class of bug one level up.

### 3.9 The laptop's brms rebuilds a cluster fit's smooths, wrongly (2026-10-02)

**Rule: post-process a fit only with the brms that fitted it, i.e. on the
cluster via `./hpc predict`, or use one of the local workarounds below.** Every
fit so far was made with brms 2.21.0, and the laptop has 2.23.1. Under 2.23.1,
`log_lik()`, `loo()`, `waic()`, `posterior_epred/linpred/predict()`,
`fitted()`, `predict()`, `conditional_effects()`, `pp_check()` and
`modelbased::estimate_*()` all return wrong numbers for these fits. **Nothing
errors and nothing warns.**

`gam_lnr6` (MullerLyer), summed over all 323,981 rows at one draw:

| draw | `lp__` | log-lik behind `lp__` | `log_lik()`, brms 2.23.1 | `log_lik()`, 2.21.0 or patched 2.23.1 |
| --- | --- | --- | --- | --- |
| 1000 (chain 2) | 38,830 | 61,083 | **-198,679** | 61,083 |
| 3500 (chain 7, the stuck one) | -55,604 | -41,083 | **-227,091** | -41,083 |

The "log-lik behind `lp__`" column was computed by hand: `cogmod::dcogmod_lnr()`
on the stored design matrices. Adding the std_normal `z`/`zs` terms (-22.2k and
-14.5k) gives back `lp__` exactly. The 2026-09-28 `gam_rdm` finding (README →
"Three traps": `mu` correlation 0.67, inverted error-rate curves) is the same
bug. Its diagnosis there, "same smooth design matrices, so not mgcv or the
BLAS", was wrong, for the reason given below.

**Mechanism.** Every post-processing function goes through
`prepare_predictions()`, whose first line is `x <- restructure(x)`. brms 2.23.0
added this to `restructure_v2()`:

```r
if (version < "2.23.0") {
  bframe <- brmsframe(x$formula, data = x$data)
  x$basis <- frame_basis(bframe, data = x$data)
}
```

This throws away the bases stored in the fit and rebuilds them from `m$data`
with the local mgcv. brms has stored those bases since 2.20.0 for exactly this
reason, so that predicting on another machine stays correct (brms#1465). 2.23.0
changed the basis format: `levels` became `group_levels`, and the GP part
changed for brms#1739. The rebuild is how old fits get the new format.

A `t2()` basis is not unique. Each `cr` marginal's penalty has a
two-dimensional null space (eigenvalues ~1e-14), and the sign of its basis
vectors is whatever LAPACK returns. Comparing the rebuilt bases with the stored
ones on the laptop gives the same pattern in every dpar (`-` is an exact sign
flip):

| matrix | columns, rebuilt vs stored |
| --- | --- |
| `Xs_<dpar>` (3) | `= - -` |
| `Zs_<dpar>_1_1` (9) | all `=` |
| `Zs_<dpar>_1_2` (6) | `= = = - - -` |
| `Zs_<dpar>_1_3` (6) | all `=` |

The draws are coefficients of the stored basis, so the flipped columns get
their coefficients with the wrong sign. That breaks two of the three
unpenalised (linear) columns and one of the three penalised components of
every smooth, at every point.

`standata(m)` and `standata(m, newdata = ...)` do **not** restructure, so they
return the stored bases. That is why the 2026-09-28 comparison of design
matrices saw nothing. It is also why `lba_modes/` is correct: those scripts
compute by hand from `standata(m)`. And it means `standata()` cannot be used to
check what `posterior_*()` computes.

**Version or platform? Both, in that order.**

- **The brms version decides whether the basis is rebuilt.** brms 2.21.0 on the
  laptop uses the stored bases and is exact (table above; `gam_rdm` matched the
  cluster to 1e-15 on 2026-09-28). Only brms >= 2.23.0 reading a fit made before
  2.23.0 rebuilds. A fit made with >= 2.23.0 is never rebuilt, but that does not
  help with the fits we have.
- **The platform decides what the rebuild gives.** Rebuilding with mgcv alone
  (`smoothCon()` + `smooth2random()`, as brms calls them) on the laptop gives
  the same flips under mgcv 1.9-0, the cluster's version (built from source
  into a scratch library), as under the laptop's 1.9.4. So mgcv's version is
  not the cause; LAPACK/BLAS is. The laptop runs R's reference LAPACK 3.12.1,
  the cluster FlexiBLAS/OpenBLAS with LAPACK 3.11.0. A rebuild on the cluster
  itself would probably reproduce the stored bases. That is untested, and
  nothing should depend on it.

**Restoring the whole stored basis is not enough either.** 2.23 reads group
levels from `basis$group_levels`, but a 2.21 basis calls them `levels`. With
`m$basis <- <stored basis>` the smooths come out right, but a `newdata` that
holds only some participants is matched to them by position. The 4
participants in a 6-row `newdata` got the `r_` of the first 4 levels in
alphabetical order. With every participant in `newdata`, or with
`re_formula = NA`, the result is correct. So the patch swaps back only the
smooths:

```r
b <- m$basis
m <- brms::restructure(m)   # what 2.23 does anyway; it also records the version, so later calls skip it
m$basis$dpars <- b$dpars    # the stored smooths; keep the rebuilt group_levels
```

`keep_stored_basis.R` wraps this as `keep_stored_basis(m)`, with
guards for the cases it does not cover. Apply it to a fit fresh from
`readRDS()`: once 2.23 has restructured a fit, the stored bases are gone.

Verified 2026-10-02 against brms 2.21.0 on the same laptop. For `gam_lnr6`,
`gam_rdm` and `gam_ddm4` (MullerLyer) the patched 2.23.1 is identical, with
max |diff| 0, on every check: each dpar's linear predictor on 500 training rows
with participants, the same on a 40-point grid without participants,
`log_lik()` on those rows, and `estimate_relation(predict = "mu")`. The
unpatched 2.23.1 on `gam_lnr6` is off by up to 9.4 on the link scale
(`sigmabias`), 5.2 (`nuone`) and 2.6 (`mu`); only `poutlier`, which has no
smooth, is unaffected. Its `estimate_relation(mu)` correlates 0.65 with the
right one, and the summed `log_lik()` over the 500 rows is -277 against +110.

The scripts behind these numbers were one-off and are not kept; this section
is the record. The patch covers what these fits contain (`t2()` smooths and one
grouping factor), not `gp()` terms, and only for fits with
`m$version$brms < 2.23.0`.

**What is affected** (audited 2026-10-02):

| where | affected? |
| --- | --- |
| `2_modelcomparison.qmd` | **no.** It reads prediction files only. All 11 in `analysis/models/predictions/` are format 2 and were made on the cluster with brms 2.21.0 (their `$meta$packages`). The local fallback, `run_predictions()`, refuses a brms other than the fit's. |
| `combine_model.R` (`loo`) | **no.** It runs on the cluster under 2.21.0. |
| `3_analysis.qmd` | **yes:** `estimate_relation()` at ~116 (the density plots) and ~307 (the heatmaps), on whatever fit it loads. It currently reads `models/gam_lnr_muller.rds`, which no longer exists, and nothing has been rendered from it, so no wrong output exists yet. Before reviving it, move its predictions into `predictions.R`, or patch the fit right after `readRDS()`. |
| `server/lba_modes/` | **no.** Computed by hand from `standata(m)` and validated against `lp__` (`lba_decompose.R`). |
| `0_preprocessing.qmd` | **no.** Its `estimate_relation()` calls are on `mgcv::bam` models fitted in the same session. |
| summaries of the draws (`summary()`, `fixef()`, `as_draws_*()`, Rhat, `lp__`) | **no.** No basis is involved. |
| ad hoc calls on the laptop | **yes**, any of the functions listed at the top of this section, plus `conditional_smooths()`, `bayes_R2()`, `insight::get_predicted()` and `marginaleffects`, on a cluster fit under 2.23.1 |

**Safe workflow, in order of preference:**

1. **Predict on the cluster**, with `./hpc predict <model>`. That uses the fit's
   own brms, mgcv and LAPACK. Anything new that needs a fit goes into
   `predictions.R` (AGENTS.md ground rule).
2. **Locally, brms 2.21.0 in a library of its own**, put first on the path
   before brms loads: `.libPaths(c("<lib>", .libPaths()))`. Build it once with
   `install.packages("https://cloud.r-project.org/src/contrib/Archive/brms/brms_2.21.0.tar.gz", repos = NULL, type = "source", lib = "<lib>")`.
   brms is pure R, so this needs no compiler. This is what the
   `run_predictions()` guard asks for.
3. **Locally under 2.23.x, `keep_stored_basis()`**, for exploratory work.
   `run_predictions()` still refuses this, on purpose: the guard compares
   versions and cannot see whether the basis was patched.
4. **Do not upgrade brms on the cluster** while these fits are in use. Pin 2.21.0
   in `install_pkgs.R`, which has been open since 3.8. A rebuild there would
   probably agree with the stored bases, but that is the same untested
   assumption that failed on the laptop.

**A check that catches it, whichever route:** at one draw, compare
`posterior_linpred(m, dpar = <any smoothed dpar>)` on a few training rows (no
`newdata`) with the same predictor computed from `standata(m)` and the draws.
Or sum `log_lik()` over all rows at one draw, add the std_normal `z`/`zs` terms,
and compare with `lp__` (`lba_modes/lba_decompose.R` does that bookkeeping).

Reported upstream on 2026-10-02, with the fix of keeping the stored `$sm`
entries in `restructure_v2()`. It reintroduces brms#1465 for every fit older
than 2.23.0, and it is still in brms master (2.23.2) as of that date.

---

## 4. Measurements

### 4.1 Dataset

The 2026-10-12 data (block exclusion on facilitating-trial errors), which the
2026-10-13 refit uses; the data of every earlier fit in brackets. The
measurements in this file were made on the earlier data, which is why they
quote 323,981 rows.

| | value |
| --- | --- |
| total rows | 962,180 (963,885) |
| total participants | 2,220 (2,221) |
| MullerLyer rows / participants | 324,580 / 2,217 (323,981 / 2,215) |
| VerticalHorizontal rows / participants (`*_verticalhorizontal` models) | 321,113 / 2,202 (319,760 / 2,200) |
| Ebbinghaus rows / participants | 316,487 / 2,181 (320,144 / 2,204) |

`fit_model.R` defaults to **all** participants since 2026-09-18 (it used to
default to 30, which risked a production submission silently fitting the test
subset). Override with `IGC_NPARTICIPANTS=30` for a smoke test.

### 4.2 Partition quotas → hub

Per-partition caps, the account-wide 550-CPU `normal` cap and the
`sussexneuro` group pool are in [hub `artemis.md#quotas`](https://github.com/RealityBending/Lab/blob/main/hpc/artemis.md#quotas).
What they mean for this project:

**Concurrency on `long` is `floor(140 / cpus-per-task)`** — with
`--cpus-per-task=16` that is **8 tasks**, whatever `--array` says, shared
across everything the account has queued there (FakeArt included: same
user). So the array widths of the jobs running at once should add up to 8.
Production is two models x `--array=1-4`; one model at `--array=1-12` would
waste a third of the throughput running a half-empty second wave.

#### 4.2.1 `sussexneuro` for this project (2026-09-22)

Access verified by bare submission (job 11404314). With the two production
arrays holding 128 CPUs at `long`'s ceiling, a 16-CPU `sussexneuro` job
(11404316) went straight to RUNNING on `artemis-rtx-02`: a further model can be
fitted alongside the `long` arrays at no cost to them, which is what
`gam_ddm5` does as of 2026-09-22 (see 6):

```bash
./hpc fit gam_ddm5_mullerlyer --partition=sussexneuro
```

Its 60-day wall is 8 usable days for us (Artemis's checkpoint request; Stan
cannot checkpoint). It does *not* on its own make full-data `gam_ddm7`
(~a fortnight per chain, 4.7.1) a good idea; the 200-participant subsample is
still the way to look at that model, and `sussexneuro` is a good home for it
because 2.7-4.7 days per chain sits inside 8 days with room to spare.

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
  and leapfrog counts, never by `sacct` elapsed alone. At full data the gap
  is ~2-4x, and every slow production shard so far was on a non-zen5 node:
  see `node_history.md` for the per-shard log, and add to it.
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
trial per rejection, and the fix is in cogmod's `NEWS.md` (0.3.3): the `A == 0`
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
divergence count. So far 8 of 8 chains at 30 participants (2026-09-18 smoke,
`gam_lnr` + `gam_ddm4`) and 6 of 6 across the three new families (4.7) started
cleanly, against a measured 7% rejection rate at that size under 0.3.2.

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

### 4.7 The new families, smoke-tested (2026-09-18)

First fits of `cogmod_rdm`, `cogmod_lba2` and the full seven-parameter
`cogmod_ddm`, all at 30 participants / 3,806 rows, warmup 300 + 100 draws,
2 chains x 4 threads on `short`. Same settings for all three, so the columns
are comparable; `gam_lnr` and `gam_ddm4` from the same matrix are included for
scale.

| model | wall | n_leapfrog | stepsize | divergent | max Rhat | params |
| --- | --- | --- | --- | --- | --- | --- |
| `gam_rdm` | 10.9 min | 129 | 0.041 | **0** | 1.119 | 276 |
| `gam_lba` | 9.5 min | 233 | 0.027 | **6%** | 1.199 | 396 |
| `gam_lnr` | 6.0 min | 267 | 0.015 | 1% | 1.321 | 336 |
| `gam_ddm4` | 7.1 min | 128 | 0.033 | 0% | 1.440 | 278 |
| `gam_ddm7` | **did not reach iteration 100 in 50 min** | — | — | — | — | — |

All three new families compile, and `cogmod_priors()` / `cogmod_inits()` handle
the hand-fixed dpars — including `gam_lba`'s `sigmazero = 1`, the only non-zero
constant in the registry. No chain was rejected at initialisation.

**`gam_rdm` is the best-behaved model in the registry**: no divergences, the
largest step size, the lowest Rhat. **`gam_lba`'s 6% divergences** are an order
of magnitude above the others here; warmup 300 is short enough that 4.6 would
predict this to clear at warmup 1000, but it is the number to check first in
its production `REPORT`.

#### `gam_ddm7` is not viable as specified

It did not reach iteration 100 in **50 minutes at 30 participants**, where the
other four finished all 400 in 6-11. CmdStan only refreshes every 100
iterations so the exact count is unknown, but netting off compile time that is
**upwards of 25x the per-iteration cost of its siblings, on 1/74th of the
production data**. Scaled against the 12-20 h per chain that 4.4.1 projects for
full data, `gam_ddm7` lands near a fortnight per chain — past `long`'s 8-day
ceiling, let alone the 2-day `--time`.

> **Corrected 2026-09-18, later the same day.** This section originally read
> "a geometry problem, not a compute problem", and guessed that the three
> between-trial variabilities were hard to identify because each carries 25
> tensor coefficients. **That diagnosis was wrong, and it is the expensive kind
> of wrong**: it makes prior tightening and reparameterising the smooths sound
> like they might work, and neither does. The cost is a *branch in the
> likelihood* and the smooths have nothing to do with it. See 4.7.1. The
> paragraph below is kept because the decision it reached still stands.

This is a geometry problem, not a compute problem, so more hours will not fix
it. The suspects are the three between-trial variabilities: `sigmadrift`,
`sigmabias` and `sigmandt` are identified through the *shape* of the RT
distribution rather than its location, weakly so even with flat predictors, and
here each carries 25 tensor coefficients plus a participant intercept. The
confirming diagnostic is `n_leapfrog` pinned at 1023 (treedepth 10), which
would say the sampler is taking maximum-length trajectories and still not
moving.

**Decision 2026-09-18: `gam_ddm7` is out of the production set**, pending a
separate investigation into how to parameterise it. `gam_ddm4` and `gam_ddm5`
are unaffected, and so are the other new families. Leave the entry in
`models.R` — it is the thing that investigation will start from.

#### 4.7.1 Why it is slow: a branch into numerical quadrature, not geometry

Measured 2026-09-18 with a standalone Stan program calling `wiener_lpdf` in
each of its forms, cost normalised by `n_leapfrog__` so the variants are
comparable. cogmod's `NEWS.md` (0.3.3) records the mechanism. The full report,
with the proposed fix and its design constraints, was
`cogmod_ddm_cost_issue.md`; it was removed on 2026-10-02 and is still in git:
`git show 402f9ae:analysis/server/cogmod_ddm_cost_issue.md`. The harness that
produced the table below was never committed and is not recoverable; the
report's §5 specifies one that would re-derive it, unrun.

`cogmod_ddm_decision_lpdf()` routes on an **exact-zero test**:

```stan
if (sw == 0 && sigmandt == 0) { ... analytic fast paths ... }
...
return wiener_lpdf(y | boundary, tau0, w, v, sigmadrift, sw, sigmandt);
```

The 7-parameter form falls back to adaptive numerical quadrature — cogmod's own
comment above that line says so. `sigmadrift` is handled analytically, which is
why `gam_ddm5` is cheap; the other two are not.

| density variant | µs/obs/gradient | x the DDM-5 path |
| --- | --- | --- |
| classic 4-parameter (`gam_ddm4`) | 1.97 | 0.36 |
| 5-parameter, `sigmadrift` free (`gam_ddm5`) | 5.56 | 1.0 |
| one of `sigmabias`/`sigmandt` free (1-D quadrature) | 100.0 | **18.0** |
| both free (2-D quadrature, `gam_ddm7`) | 307.1 | **55.3** |
| both free but driven to `1e-5` | 144.7 | 26.0 |

The benchmark reproduces the 4.7 table it is explaining: it predicts
`gam_ddm4` at 5.2 min against the 7.1 measured, and `gam_ddm7` reaching ~9
iterations in 46 min against the observed "did not reach 100". At full data
with 8 threads the classic path predicts 80 ms per gradient against the 78 ms
4.4.1 measured for `gam_lnr`.

**What this rules out, and what it leaves open:**

- **Tight priors near zero do not work.** 55x → 26x and no further: `1e-3` and
  `1e-5` measure the same, because the fast path tests for *exact* zero and an
  estimated parameter never satisfies it. Note also that "zero" is on the link
  scale — `sigmabias` is logit (η = 0 is *half* the maximum start-point range)
  and `sigmandt` is log (η = 0 is **one second** of non-decision-time range).
  A prior centred at 0 there is pathological, not conservative. cogmod's own
  `normal(-2, 1)` and `normal(-3, 1)` are already the sensible locations.
- **The smooths are innocent.** The cost is per-observation and independent of
  how many coefficients feed the linear predictor, so making the three
  variabilities intercept-only changes nothing. Do not spend a run finding
  this out.
- **A warm start does not close the gap either.** It buys back warmup, not
  per-gradient cost; 500 retained draws alone are ~18 days at full data.
- **Freeing exactly one variability costs 18x, not 55x** — the same either way
  (18.0 against 18.0 to three figures), so which one is a modelling choice,
  not a performance one. If one is ever freed it should be **`sigmabias`**:
  freeing `sigmandt` silently redefines `ndt` as the *lower bound* of a
  between-trial distribution rather than its midpoint (cogmod documents this),
  and `ndt` carries its own 2-D smooth in every model in the registry, so the
  `ndt` surface would stop being comparable across the model set.
- **Subsampling works, and is the only thing here that does.** The cost is per
  row, so per chain at warmup 1000 + 500 draws, 16 threads, 1 chain per task:

  | | full data | fits ~5 days on `long` |
  | --- | --- | --- |
  | `gam_ddm7` as specified | ~54 days | ~240 participants |
  | `gam_ddm7`, the two pinned tight | ~26 days | ~500 participants |
  | a six-parameter DDM | ~18 days | ~700 participants |

  At 200 participants `gam_ddm7` is 2.7 days per chain if leapfrog settles at
  255 and 4.7 at 511, against 8.7 — over the wall — if treedepth stays pinned
  at 10. That last case is itself the finding, so the run is informative
  whether or not it completes. Give it **1 chain x 16 threads** (4.4.1:
  threading scales, more chains per task does not) and **its own
  `IGC_MODELS_DIR`**.
- **The real fix is in cogmod**, not here: replace the delegation to Stan's
  adaptive quadrature with cogmod's own fixed-node Gauss-Legendre rule, the
  pattern already used for the inverse-Gaussian's drift integral
  (`.pwald_sv()` / `.WALD_STAN_PRELUDE` in `R/core_shifted.R`). A 3x3 rule
  projects to ~6x faster, which is what would make a six-parameter DDM routine
  and a seven-parameter one reachable in combination with a warm start.
  The design constraints are §4.1 of the removed report (the `git show`
  above).

### 4.8 Criterion cost — `loo` over every draw is affordable (2026-09-20)

Measured on the 3-shard `gam_lnr` combined fit (6 chains x 500 = **3,000
draws** over **323,981** responses) on `general`, 16 CPUs, 256 GB requested.
Elapsed seconds for `brms::waic()` / `brms::loo()`, the calls `add_criterion()`
makes:

| draws | `waic` | `loo` |
| --- | --- | --- |
| 500 | 169 s | 300 s |
| 1500 | 229 s | 435 s |
| 3000 (all) | 487 s | **569 s** |
| 500, `cores = 16` | — | 328 s |
| 3000, `cores = 16` | — | 535 s |

Reading and merging the three shards is 8 s of that. The log-likelihood matrix
is 3.6 GB at 1500 draws and 7.2 GB at 3000.

**Memory: 320G since 2026-10-09**, because fits now carry 8,000 draws (6). The
peaks below scale linearly with draws, which gives ~210 GiB for the
5-6-parameter models and ~245 GiB for `gam_lnr6` / `gam_lba`. This is a
projection: measure the first 8,000-draw combine. A 4,000-draw combine still
fits in `--mem=192G`.

**Memory: 192G since 2026-09-29, not 128G.** The 128G was set from the 45 GB of
a `waic`-at-1,500-draws combine. With `loo` over every draw of four shards, the
peak is several times that. It scales with draws x rows: `log_lik()` evaluates
each distributional parameter over all draws x rows, and `loo()` makes several
copies of the log-likelihood matrix, each 9.7 GiB at 4,000 draws x 323,981 rows.
`MaxRSS` (cgroup, sampled every 30 s unless noted):

| combine | shards (draws) | peak |
| --- | --- | --- |
| `gam_ddm4`, `gam_lnr6`, `gam_ddm5` | 2 (2,000) | 53-65 GiB |
| `gam_lnr`, `gam_ddm4`, `gam_ddm5`, `gam_lnr_vh`, `gam_ddm4_vh` | 4 (4,000) | 100-106 GiB |
| `gam_lnr6` (11417254, 128G) | 4 (4,000) | OOM-killed at 17 min, in `log_lik()` |
| `gam_lnr6` (11417275, 192G, 5 s sampling) | 4 (4,000) | **122 GiB**, 18.5 min |

`gam_lnr6` has seven draws x rows parameters (one more than `gam_lnr`), and at
122 GiB resident it sits at 95% of 128 GiB before page cache and the spikes
between samples. Expect the same from `gam_lnr6_*` on the other illusions and
from `gam_lba*`, which also have seven. The 30 s rows understate their true
peaks too, so 128G was never a comfortable margin at four shards. The general
nodes had 300-990 GB free, and the 192G job dispatched within two minutes.
Ignore the 185 GB peak of the probe job: it ran all five cells in one R
session and kept every result.

Three measurements set the defaults:

- **Nothing is near the wall.** The dearest cell is under 10 minutes against
  `general`'s 8 h (4.2). `combine.slurm`'s 128 GB and its original 6 h were
  reasoned from an assumption that the custom family's R-level `log_lik` would
  dominate a serial loop over 324k responses. It does not.
- **Cost is strongly sublinear in draws.** `loo` over all 3,000 draws costs
  **1.9x** what it costs over 500 — for **6x** the draws. Most of the price is
  fixed overhead, not per-draw work, so subsampling buys very little wall clock
  while throwing away precision.
- **`cores` does not help, and can hurt.** It does reach `log_lik()` —
  `add_criterion()` passes `...` to `waic.brmsfit`/`loo.brmsfit` and on to
  `log_lik.brmsfit`, which takes it — but at 3,000 draws `cores = 16` moved
  569 s to 535 s (6% for 16x the CPUs), and at 500 draws it moved 300 s to
  **328 s**, i.e. 9% *slower*: the fork and the copies cost more than the work
  they split. Do not reach for it again, and do not size `combine.slurm` on
  the expectation that it will pay — 4 CPUs is deliberate.

Hence `combine_model.R` defaults to `loo` over every draw. `loo` is also the
criterion the comparison wants on its own merits: PSIS-LOO carries per-observation
Pareto-k diagnostics that say when the estimate is untrustworthy, and WAIC has
none. `waic` is 1.2-1.5x cheaper at equal draws and stays available through
`IGC_CRITERION=waic`, as does a draw count through `IGC_CRITERION_NDRAWS` — for
a model that actually overruns, not pre-emptively.

#### Pass no `ndraws` rather than `ndraws = ndraws(m)`

They are not the same call. Omitting the argument leaves the draws in their
chains, which is what `loo`'s `r_eff` needs in order to discount
autocorrelation; passing a count routes through the subsampling path, which
flattens the chain structure and degrades `r_eff`. Since every draw is the
default, `combine_model.R` omits the argument on that path rather than passing
the total.

#### Shards are kept by default now, and that reopens 3.6

Deleting shards after a successful combine was a guard, not housekeeping:
`file_refit = "never"` means a shard left on disk is silently reused even when
the formula or the data changed (3.6), and deleting kept "the file exists"
equivalent to "this fit is finished and banked".

It is still worth keeping them. A combined fit is ~10 minutes to rebuild from
shards — a different criterion, a re-run after a cogmod change — while each
shard is ~43 h of compute (4.4), against disk that is free. The 2026-09-20
combine would otherwise have destroyed three shards that were about to be
re-combined under `loo`.

So the guard is now manual, and it is the one thing to remember about this
change: **a changed parametrisation or formula needs its own `IGC_MODELS_DIR`,
or `IGC_FILE_REFIT=always`.** `IGC_DELETE_SHARDS=1` restores the old behaviour.

### 4.9 `gam_lba` (MullerLyer) did not converge: two modes (2026-10-01)

Fitted and combined on the second account (4 shards x 2 chains, warmup 1000 +
500, brms 2.21.0, cold starts, the post-5.7 priors with `exponential(1)` on
`sds`). The combined file was uploaded to `models/combined/` on `dmm56`. **It
has no prediction file and is not in `2_modelcomparison.qmd`.** Do not run
`./hpc predict gam_lba` on it as it stands.

The eight chains sit in two modes and none of them crosses between them. lp__
is flat from the first 100 draws to the last in every chain:

| | chains (shard.chain) | lp__ | lprior | `driftone` | `sigmaone` | `sigmabias` | `boundary` | sd(`driftone` \| ppt) | sd(`boundary` \| ppt) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | 1.1, 1.2, 2.2, 3.2, 4.2 | ~36,965 | ~ -91 | -4.9 | 2.0 | -0.04 | 0.38 | 1.23 | 0.55 |
| B | 2.1, 3.1, 4.1 | 38,226-38,308 | ~ -143 | -7.0 | 2.6 | -0.45 to -0.78 | 0.57-0.64 | 1.77-1.86 | 0.34 |

The parameter columns are intercepts on the link scale. `driftone` is identity
and the other three are softplus. The `mu`, `ndt` and `poutlier` intercepts
barely differ between the modes.

The investigation is written up in **`cogmod_lba_modes_issue.md`**, and the
scripts are in `lba_modes/`. In short:

- **The modes are real.** Along the straight line between them the
  log-likelihood dips by 4,500-5,700. A is a tight local optimum: its 5 chains
  agree (Rhat 1.03). B fits better by **1,364 log-likelihood units** (lp__ split
  by hand, then checked observation by observation), but its 3 chains disagree
  (Rhat 1.86, on `ndt` / `sigmabias` / `boundary`).
- **B's gain is in correct-response RTs** (+1,433, against -69 on errors),
  mostly in conflicting D2-D4. It is spread across participants, not
  concentrated in a few.
- **cogmod's Stan density is accurate** (to 1e-11 against direct integration,
  with smooth gradients), so the code is not the cause. The geometry is.
  Facilitating trials have 1-3% errors, so the error accumulator sits at
  v/s of -4 to -5.4. There the truncated drift becomes Exponential-like and is
  identified only through |v|/s². The modes place it at drift -7 to -12 and
  -15 to -17 with the same realised drift, and they split threshold against
  start-point range differently across the design.
- **Both modes predict the same behaviour**, and both miss fast errors in the
  easy cells the same way. A longer run would not reveal the "true"
  parameters.
- **The diagnostics show the split:** max Rhat 13.8 on a `sigmabias` smooth
  coefficient, and population SDs at 7-10, against the 1.3-1.9 that
  participant terms give the other models (6.1). 6,342 of 15,690 parameters
  have Rhat > 1.1. There were no divergences, and 492 observations have Pareto
  k > 0.7, which is ordinary here. Divergences and k do not show this problem.
- **Every iteration of every chain hits treedepth 10** (1,023 leapfrogs). Step
  sizes are 0.0017-0.0052. Chains took 48-72 h each (31-50 h of warmup), so a
  shard runs 2.3-3 days.
- **The priors are overwhelmed, not causal.** The `driftone` intercept sits
  3-4 SDs below `normal(1, 2)`, and the smooth SDs are 15-36 against
  `exponential(1)`. The prior difference between the modes (52) is small next
  to the likelihood difference (1,364).
- **The init is worth fixing in cogmod:** LBA2 starts `driftone = mu = 3`, the
  cold-start pattern cogmod already corrected for the RDM. But all 8 chains
  started at the same point and still split, so it is unlikely to be the
  whole story.

More shards of the same model would only add more chains to A or B. A warm
start near mode B would also do, but 5.5 rules it out.

**Plan (user's decision, 2026-10-01):**
1. **First, adjust cogmod's defaults**, listed in
   **`cogmod_lba_priors_inits.md`**: the LBA2 `driftone` start, the `sds` start
   and prior on the drift smooths, and the `mu`/`driftone` prior asymmetry.
   Then the user re-runs `gam_lba` unchanged to see whether the split
   replicates. Give the re-run its own `IGC_MODELS_DIR` (or
   `IGC_FILE_REFIT=always`), and `./hpc install cogmod` first.
   *Status 2026-10-02:* this is done in cogmod 0.3.4 on `dev` (2), apart from
   the warning that B3 of the note proposes; the note's header lists what
   was changed. The next step is the re-run itself, on 0.3.4 like every new
   fit (2).
   *Status 2026-10-09:* **re-run done, and the split replicated.** The
   `gam_lba_MullerLyer.rds` on the laptop is a 0.3.4 fit (`732be30`, sampled
   2026-10-07). It is a new run, not the old one re-saved: the draws differ
   from `gam_lba_MullerLyer_20261001.rds`. It shows the same 5/3 clusters at
   the same lp__ levels (about 37,000 against 38,300), and the
   VerticalHorizontal fit splits too (4.10). That triggers step 2.
2. **Only if it replicates**, change the model:
   **`sigmaone ~ 1 + (1 | Participant)`** (no smooth on `sigmaone`). With the
   error accumulator's drift SD shared across conditions, |v|/s² pins
   `driftone` in every condition. It would be a new `models.R` entry, not an
   edit to `gam_lba`. Further options, in reserve: `sigmabias` treated the same
   way, `sigmaone = 1`, or dropping the LBA in favour of `gam_lnr6`
   (`cogmod_lba_modes_issue.md` §6). Note that `gam_lnr6` has the same
   problem (4.10).

### 4.10 Convergence across the registry: three kinds of failure (2026-10-09)

Every combined fit on the laptop was examined chain by chain, from the draws
only (safe with the laptop's brms, 3.9). That covers MullerLyer and
VerticalHorizontal for all families. Ebbinghaus and `gam_rdm5_verticalhorizontal`
were seen only through the summaries stored in their prediction files. The
checks were: per-chain lp__ by quarter of sampling; Rhat on the first half of
the draws against the second; and the distance between chains, as the RMS over
the population terms of the gap between chain means, in within-chain SDs.
Sampling noise alone gives 0.1-0.3.

**No fit is still moving.** In every fit, each chain's lp__ is flat from the
first quarter of its draws to the last, and Rhat on the first half equals Rhat
on the second. So no group is "almost converged": a longer warmup or more
draws would not change the picture in any of them. What separates the fits is
whether the chains agree.

#### A. Converged at the population level: DDM-4, DDM-5, LNR, RDM

All 8 chains agree (distance 0.1-0.3). The population Rhat is at most 1.04
(1.10 for DDM-4 VerticalHorizontal) and the lp__ Rhat at most 1.02. The
Ebbinghaus DDM-4, DDM-5 and RDM summaries look the same. What remains:

- **Low ESS on a few population terms**, mostly `b_ndt_Intercept` (bulk ESS
  53-310 over 4,000 draws). More draws would fix this; a longer warmup would
  not.
- **3-25 participants per fit (66 in all) that do not mix**, with Rhat 1.1-1.8 on
  their own terms. Each such participant's chains sit in one of two tight
  solutions and never switch. Examples: DDM-4 VerticalHorizontal, Validation
  S42, `ndt` -0.035 in 4 chains and +0.089 in the other 4 (within-chain SD
  0.01). DDM-5 MullerLyer, EEG sub-4, `poutlier` 5.0 against -0.5, i.e.
  whether that participant's trials are treated as contaminants. The parameter
  is always `ndt`, usually traded against `boundary` and/or `poutlier`. Ten
  participants recur across models. RDM has none.
  `mixing_participants.R` writes the list (`poorly_mixing_participants.csv`),
  and `2_modelcomparison.qmd` → Convergence → "Participants That Do Not Mix"
  shows them.

  **Their profile.** Compared participant by participant with the others, they
  are faster (median RT at the 22nd percentile of the others), more consistent
  (RT IQR at the 15th), more error-prone (error rate at the 85th), and their
  errors are no slower than their correct responses (22nd).

  **There is no rule for excluding them a priori.** Every threshold on those
  features removes far more participants who mix well than ones who don't.
  For example, "IQR < 0.12 s and errors > 30%" catches 8 of the 71 flagged
  participant-illusion pairs and 74 others. "Any RT < 0.15 s" catches 8 and
  201. It is the tail of a continuum, not a separate group.

  **Decision (user, 2026-10-09):** no participant is flagged or excluded
  beforehand. The models are refitted with longer runs. Any participants who
  still do not mix are reported as a result in their own right: a response
  profile these models cannot accommodate. The evidence above predicts that
  the longer run will not merge the per-participant solutions. It does give
  each chain more chances to switch, and the refits on 0.3.4 are due anyway
  (6.4). Rerun `mixing_participants.R` on the refits to get the final list.

#### B. Separate stationary modes: every model that frees `sigmabias`

LNR-6, RDM-5 and LBA, on both illusions. Their worst-mixing parameters are
`sigmabias` terms (intercept, smooth, `sds`, participant SD and participant
terms), along with whatever trades off against `sigmabias`: `mu`'s smooth in
the LNR and `driftone` in the LBA. No model with `sigmabias = 0` has this
problem. The chains form clusters. Within a cluster they agree (distance 0.1);
between clusters they are 3-6 SDs apart on average over all the population
terms. No chain crosses between clusters:

| fit | clusters (combined chain index) | lp__ of the other clusters, against the largest |
| --- | --- | --- |
| `gam_lba` ML, 0.3.4 re-run | {1,4,6,7,8} / {2,3,5} | +1,350: the 4.9 split replicated |
| `gam_lba` VH, 0.3.4 | {3,5,6,7} / {2,4,8} / {1} | +1,560 / -70 |
| `gam_rdm5` ML | {1,2,3,4,5,8} / {6} / {7} | **+1,670** / -70 |
| `gam_rdm5` VH (summary only) | 5 chains / 3 chains | -2,300 |
| `gam_lnr6` VH | {1,6,7,8} / {2,3,4,5} | about -80; participant `sigmabias` changes sign between modes |
| `gam_lnr6` ML | {2,4,5,6,8} / {1} / {3} (+ chain 7 frozen, C) | -1,015 / about 0, 3 SDs away |

Points worth keeping:

- **The modes are local traps, not one wide posterior.** They differ by
  hundreds to thousands of log-density units, and the cluster with the most
  chains is often not the best one: in `gam_rdm5` ML a single chain out of 8
  found a mode 1,670 units higher. There may be modes no chain has found.
- `gam_lnr6` ML chain 3 is the clearest case of a ridge. It matches the main
  cluster's lp__ but sits 3 SDs away.
- Most of these models hit treedepth 10 on every iteration (1,023 leapfrogs);
  `gam_lnr6` VerticalHorizontal does not. Over 500 draws no chain moved. More
  warmup or more draws will not join the clusters.
- **Their prediction files mix modes.** `gam_lnr6` and `gam_rdm5` are plotted
  in `2_modelcomparison.qmd`, and their curves and `loo` average over every
  cluster. Treat them as unreliable until refitted.
- **Exception: `gam_rdm5_ebbinghaus` (cogmod 0.3.4, combined 2026-10-13)
  converged.** All 8 chains sit within 20 lp__ units of each other (SD ~115),
  with population Rhat 1.010, lp__ Rhat 1.007, minimum bulk ESS 761 and 0.02%
  divergent. Its `loo` is ~3,000 elpd above RDM's on the same illusion. So the
  split is not universal to a free `sigmabias`. Whether that comes from 0.3.4
  (the ML and VH `gam_rdm5` fits are 0.3.3) or from the illusion is open:
  `gam_lnr6_ebbinghaus` (0.3.4) does not answer it: it failed earlier, in
  adaptation, like `gam_lnr_ebbinghaus` (C below).
- **This is a result in itself**, to report in the paper. On this data, with
  trial-level smooths on every parameter, the start-point range cannot be
  identified by any racing family. Freeing it turns a well-behaved model
  (RDM: 0 participants above Rhat 1.05) into a multimodal one.

**Decision (user, 2026-10-13): first refit the full formula as it is**, with
every other model, under the new harness (cogmod 0.3.5, the 2026-10-12 data,
1500 warmup + 1000 draws, 4.11). Only if the split replicates there is a
simplified variant fitted, as a new `models.R` entry, and the failure of the
full formula is then reported in the paper as a finding in its own right
(the Discussion points below). The options for that case (each a new
`models.R` entry, never an edit to an existing model):

1. **A simpler formula: `sigmabias ~ 1 + (1 | Participant)`**, so no smooth on
   the start-point range. Its smooth and participant terms are where the
   clusters disagree most (`sds_sigmabias` at 2, 12 and 22 in the three
   `gam_rdm5` ML clusters). This is the same move already planned for the
   LBA's `sigmaone` (4.9). A further step is `sigmabias ~ 1`, a single range
   for everyone. **This is the preferred option** (2026-10-09), for two
   reasons the paper's Discussion is to make (TODO in
   `paper/manuscript.qmd`):
   - *Theory: a motivation, not a proof.* The start point precedes the
     stimulus, and trials are interleaved. But the DDM's `bias` is smoothed
     too and does vary with the predictors, so this argument alone cannot
     single out `sigmabias`. What separates them: with `dec(Error)`
     (accuracy coding), a pre-stimulus `bias` would average 0.5, so a
     stimulus effect on it is early stimulus-driven processing (e.g. an
     automatic pull of the illusion), and that is directional. `sigmabias` in
     the race models is a symmetric range shared by both accumulators, with
     no direction. Their directional early pull lives in `mu` vs `driftone`.
     Its DDM counterpart is start-point *variability*, fixed at 0 in DDM-4 and
     DDM-5. So an effect on `bias` does not predict one on `sigmabias`.
   - *Identifiability, the decisive argument.* The response depends on the
     distance from start point to boundary, and `boundary` already varies
     with the conditions. A condition effect on `sigmabias` is nearly
     interchangeable with one on `boundary`, which is the trade-off that
     separates the modes. The empirical contrast: DDM-4 and DDM-5 converge
     with `bias` smoothed, while every race model with `sigmabias` smoothed
     does not.

   Linear terms (`~ Illusion_StrengthZ + Illusion_DifferenceZ + (1 |
   Participant)`) would shrink that trade-off from about 24 coefficients to 2,
   but not remove it. A single slope on the *signed* strength would also force
   the two sides of the illusion to mirror each other. That variant is worth
   fitting only as an explicit test of selective influence (slopes near 0),
   not as a fix. `(1 | Participant)` is kept because people plausibly differ
   in start-point variability, but its terms also disagree between modes (sign
   flips in `gam_lnr6` VH), so `~ 1` is the fallback if it still splits.
2. **Tighter priors on `sigmabias`.** This is less promising. In the LBA the
   likelihood gap between modes is 25 times the prior gap (4.9), so only a
   prior close to fixing the parameter would choose between them. Note that
   0.3.4's `sds` prior is *looser* than 0.3.3's on these `t2()` smooths
   (5.7). The `gam_lnr6` and `gam_rdm5` fits in the table are 0.3.3, and the
   two `gam_lba` fits are 0.3.4, so cogmod's 0.3.4 defaults did not remove
   the split. The 0.3.4 Ebbinghaus runs of `gam_rdm5` and `gam_lnr6` (6.4)
   will show the same for those two. A tighter-prior variant would set the
   `sds` prior on `sigmabias` explicitly, rather than through cogmod's
   default.
3. **Test either option cheaply first.** The question is whether 8 chains
   still split on a subsample of 200-480 participants with the current
   formula. If they do, the variants can be compared at that size before any
   full-data run.

#### C. Frozen chains

A chain whose step size collapsed during adaptation. Its lp__ barely moves, far
below the others, and its smooth coefficients stay pinned near their start.
This is a sampler failure, not a mode, and a restart of that shard fixes it.

- `gam_lnr6` MullerLyer, combined chain 7 (shard 4, chain 1). Step size 3e-5
  against about 8e-3, lp__ -94,600 below the other chains, lp__ SD 42 against
  about 120.
- `gam_lnr_ebbinghaus` on 0.3.3: step sizes down to 1e-8, and lp__ 40,000-82,000 below
  the others in two chains. Its 0.3.4 refit failed the same way (next item).
- **The 0.3.4 refit of `gam_lnr_ebbinghaus` (11553388) failed in all 8
  chains, and so did `gam_lnr6_ebbinghaus` (11553384)**: adapted step sizes
  1e-4 to 6e-4 (one chain 9e-6), treedepth 10 on every draw, max Rhat 3.0,
  and each chain drifting monotonically through its 500 sampling draws in the
  scale parameters (`lnr_modes/lnr_ebb_chains.R`). That is neither a mode
  split nor a single frozen chain: a warmup that never reached the typical
  set, in both LNR models on this illusion, while `gam_rdm5_ebbinghaus`
  converged on the same data. The 1500-warmup refit (4.11) is the test of
  whether more warmup is enough; if not, reparametrise rather than lengthen
  again.

Restarting `gam_lnr6` MullerLyer's shard 4 alone is not worth the compute while
the model is in group B. `combine_model.R` now prints one line per chain (step
size, treedepth, divergences, lp__ level and spread). It warns on a frozen
chain, naming the shard, and on a chain whose lp__ is more than 5 within-chain
SDs from the median, which is likely a separate mode. Chains in the same mode
differ by about 0.3 SD. A mode at the same lp__ (`gam_lnr6` ML chain 3, or
`gam_lnr6` VH) shows only in the parameters, so the check cannot see it.
Tested on `gam_lnr6` ML: it finds chain 7 as frozen and chain 1 (8.5 SDs) as
off.

### 4.11 Run length for the refit: 1500 warmup + 1000 draws (2026-10-13)

Every model is refitted on cogmod 0.3.5 and the 2026-10-12 data (blocks
excluded on facilitating-trial errors). The run length was sized from what
the finished shards measured, not from the 4.4.1 projections. Each `.out`
carries CmdStan's own split (`REPORT time per chain`: the `$chains` rows,
below the `$total NA` line that a short grep stops at), so the cost of a
warmup iteration and of a sampling iteration is known separately for every
shard. On zen5 nodes, 2 chains x 8 threads, 1000 + 500, seconds per
iteration, median [slowest chain]:

| model | warmup s/it | sampling s/it | chain wall at 1000 + 500 |
| --- | --- | --- | --- |
| LNR ML / VH / Ebb | 57 / 42 / 90 [77 / 58 / 156] | 54 / 34 / 76 [66 / 42 / 137] | 24 / 17 / 35 h |
| DDM-4 ML / VH / Ebb | 70 / 55 / 53 [87 / 70 / 115] | 83 / 39 / 35 [104 / 55 / 39] | 31 / 21 / 20 h |
| DDM-5 ML / VH / Ebb | 84 / 56 / 61 [93 / 60 / 85] | 59 / 43 / 70 [61 / 59 / 74] | 31 / 22 / 27 h |
| RDM VH / Ebb | 67 / 59 [76 / 80] | 88 / 66 [108 / 96] | 30 / 27 h |
| RDM-5 VH / Ebb | 120 / 129 [156 / 178] | 158 / 174 [181 / 223] | 56 / 61 h |
| LNR-6 ML / VH / Ebb | 180 / 137 / 202 [234 / 156 / 221] | 148 / 100 / 164 [208 / 112 / 175] | 72 / 51 / 79 h |
| LBA (colleague's account, 4.9) | ~144 [180] | ~144 [158] | 60 h |

Three things follow:

- **A warmup iteration costs 0.75-1.5x a sampling one**, so wall time is
  close to proportional to warmup + draws. The 4.4.1 picture of a short,
  expensive transient followed by cheap iterations does not hold at full
  data: most models sit at treedepth 9-10 throughout.
- **The registry splits into a fast family (LNR, DDM-4, DDM-5, RDM: 35-100
  s/it) and a slow one (LNR-6, RDM-5, LBA: 100-235 s/it)**, which are also
  the converged and the mode-splitting models (4.10).
- **Non-zen5 nodes are 2-4x slower** (a40-11 ran DDM-5 at 250-420 s/it), so
  the projections below assume `--constraint=amd_zen5`.

Projected wall per shard in days, median chain and slowest chain seen
(+0.1 day for compile, data and save):

| model | 1000+1000 median | 1000+1000 slowest | **1500+1000 median** | **1500+1000 slowest** | 2000+1000 median | 2000+1000 slowest |
| --- | --- | --- | --- | --- | --- | --- |
| LNR, DDM-4, DDM-5, RDM | 1.0-2.0 | 1.3-3.5 | 1.2-2.5 | 1.6-4.4 | 1.5-3.1 | 1.9-5.3 |
| RDM-5 | 3.3-3.6 | 4.7 | 4.0-4.4 | 5.8 | 4.7-5.1 | 6.8 |
| LNR-6 | 2.8-4.3 | 5.2 | 3.6-5.5 | 6.6 | 4.4-6.7 | 7.9 |
| LBA | 3.4 | 4.0 | 4.3 | 5.1 | 5.1 | 6.1 |

**Decision (user, 2026-10-13): warmup 1500, 1000 draws, for every model.**

- 2000 warmup puts the slowest LNR-6 and RDM-5 chains at 6.8-7.9 days against
  `long`'s 8. A shard that hits the wall is lost outright: there is no
  checkpoint, and `file_refit = "never"` only reuses finished shards. The two
  `gam_lnr6_ebbinghaus` shards still running on 2026-10-13 had taken 95 h for
  1500 iterations. If 2000 is ever wanted, run it on `sussexneuro` (60 days).
- 1500 draws instead of 1000 would need about 480G for the combine (4.8),
  which only the highmem nodes have.
- The whole refit is about 4,000 CPU-days at 1500 + 1000 (3,200 at 1000 +
  1000, 4,750 at 2000 + 1000), the nine slow fits being 63% of it: 28 days of
  wall on `long`'s 140 CPUs alone, roughly 2-3 weeks with `sussexneuro`
  alongside.

What the extra warmup can and cannot do. The converged models were stationary
from their first draw (4.10 A), so for them the gain is the draws. The
mode-split models will not be joined by any run length (4.10 B); the decision
there is to give the full `sigmabias` formula one more go under the new
harness before any simplified variant. The one failure a longer warmup
directly targets is the Ebbinghaus LNR pair (4.10 C): chains that still drift
through the whole sampling phase, i.e. a warmup 1000 that never reached the
typical set. If 1500 does not fix that, the answer is a reparametrisation,
not more hours.

The per-chain seconds come straight from the `fit_*.out` files in scratch;
one grep over their `$chains` rows reproduces the tables.

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
trade-off. brms also offers `"on_change"`, which refits when the formula,
data or prior change. It is no longer supported here, because since
2026-10-06 `fit_model.R` does the reuse and the save itself (5.3.1).

#### 5.3.1 Saving a shard without losing it (2026-10-06)

On 2026-10-05, five Ebbinghaus shards finished sampling and then died in
`saveRDS` with `error writing to connection`: `gam_rdm_ebbinghaus` 1 and 4,
`gam_ddm5_ebbinghaus` 1-3. Each had run 23-33 h. Each left a truncated `.rds`
(68K to 114M, where a good shard is 166M or 203M) and no draws: `brm(file = )`
wrote straight to the shard's final name, and the Stan CSVs were in the
node's tempdir, deleted with the job. The account was at 68 GB with no quota
set, and `/mnt/lustre` was 93% full (1.4 PB). Shards saved at 13:05 and 21:41
that same day are fine. So this was the shared filesystem failing for a few
hours, not our usage, and tidying our files cannot prevent it.

`fit_model.R` now:
- writes the Stan CSVs to `models/stan_csv/<shard>_<jobid>/` and deletes them
  only once the shard is saved. After a failed save, rebuild the fit from them
  (the recipe is in the script) instead of sampling again;
- saves to `<shard>.rds.tmp` and renames it, so a failed write leaves no file
  under the shard's name;
- reuses an existing shard only if `readRDS()` can read it, and refits an
  unreadable one, as brms itself did.

There is no automatic retry. If a save fails, the job fails and sends its
FAIL mail, and the CSVs are still there.

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
`analysis/2_modelcomparison.qmd` that writes it are untouched.

### 5.6 Everything is parameterised by environment variable

No user paths are hard-coded; `./hpc` passes them via `sbatch --export`. The
variables and their defaults are in `README.md` → "Paths" and
"Run-shaping variables" (one list, so it cannot drift).

Every fit prints one `REPORT <model> | wall | n_leapfrog | treedepth | stepsize |
divergent | accept | max Rhat | min neff_ratio` line to the `.out` log, so runs
can be compared without pulling the `.rds`.

### 5.7 The slope prior fills flats only; it does not overwrite cogmod's

`fit_model.R` adds a blanket normal on class `b` because brms leaves slopes flat
(`normal(0, 1)` until 2026-10-12, `normal(0, 2)` since; `cogmod_priors_review.md`).
Until 2026-09-18 it did that **unconditionally**, with `replace = TRUE`, which
overwrote the slope priors cogmod sets for the dpars it knows to be hard to
identify. Measured over the registry, before and after the loop:

| dpar | cogmod sets | the old loop made it |
| --- | --- | --- |
| `mu`, `bias`, `boundary` (DDM) | *flat* | `normal(0, 1)` — correct, keep |
| `ndt` | `normal(0, 0.2)` | `normal(0, 1)` — **5x wider** |
| `sigmadrift`, `sigmabias`, `sigmandt`, `sigmaone`, `boundary` (race) | `normal(0, 0.5)` | `normal(0, 1)` — **2x wider** |
| `driftone` (`gam_rdm`, `gam_rdm5`) | `normal(0, 2)` | `normal(0, 1)` — **2x narrower** |
| `driftone` (`gam_lba`) | `normal(0, 1.5)` | `normal(0, 1)` — **1.5x narrower** |

So it was not a uniform loosening: it loosened the shift and the variability
parameters and **tightened the second accumulator's drift**, which is an effect
of interest rather than a nuisance term. The loop now skips any dpar whose
blanket `b` row cogmod already filled.

**The production runs launched on 2026-09-18 predate this fix** and carry the
old priors. Whether that is worth a restart is a judgement call and not one to
make silently — the variability terms are merely less regularised than
intended, but `driftone` in `gam_rdm` / `gam_rdm5` / `gam_lba` is shrunk
harder than cogmod intends, and those are the models where it matters.

One thing the loop still does **not** touch: `sds`, the smooth wiggliness SD.
Left on brms's `student_t(3, 0, 2.5)` it is the loosest prior in the model — on
a logit or log link a half-t with median ~1.9, enough for the smooth alone to
walk a `sigmabias` across its whole range or move a `sigmandt` by a factor of
seven, undoing the tight intercept the family set on purpose.

**cogmod does set it, as of 0.3.3** — this section previously said it did not,
which was true of 0.3.2 and is not true of what the cluster now runs. The bug
was subtler than a missing prior: brms fills the *blanket* `sds` row itself and
leaves the per-term rows empty, so cogmod's filler never saw a candidate and an
`sds` on a dpar silently kept brms's default while `?cogmod_priors` advertised
`exponential(1)`. 0.3.3 replaces the blanket row instead
(`R/cogmod_priors.R:920-933`; cogmod's `NEWS.md`, 0.3.3, has the detail).

Verified 2026-09-22: the cluster library is cogmod 0.3.3 at commit `d04c7f8`,
built 2026-09-20, and that commit contains the fix. So **fits submitted from
2026-09-20 onwards carry `exponential(1)` on `sds`; anything built before that
does not** — which makes it a property of the library at submission time rather
than of anything in this repo. Worth knowing before comparing smooth surfaces
across models fitted at different times.

cogmod 0.3.4 changes this again (2). The fixed `exponential(1)` becomes a rate
set per smooth from its basis. On these `t2()` smooths one unit of `sds` is
only about 0.1 link units of wiggle, so 0.3.3's `exponential(1)` was tight
there, and 0.3.4's prior is much looser. Fits built after a 0.3.4 install are
a third prior set.

---

## 6. The production configuration

Settled 2026-09-18 from the measurements above. `fit.slurm`'s `#SBATCH`
defaults *are* this configuration, so the production run is two bare commands:

```bash
./hpc push
./hpc fit gam_lnr_mullerlyer
./hpc fit gam_lnr6_mullerlyer
./hpc fit gam_ddm4_mullerlyer
```

| setting | value | why |
| --- | --- | --- |
| partition | `long`, `--time=8-00:00:00` | `long`'s actual max (3.7) — 12-20 h per chain was the projection (4.4.1), but three `gam_lnr` shards measured 42.7-44 h on 2026-09-20, so request the ceiling rather than the guess; `--time` doesn't affect dispatch here (4.2's CPU cap does) |
| participants | `IGC_NPARTICIPANTS=all` (default) | MullerLyer: 2,217 participants, 324,580 rows since 2026-10-12 (2,215 and 323,981 before; 4.1) |
| warmup | `IGC_WARMUP=1500` (1000 until 2026-10-13) | 300 provably does not adapt at full data (4.4.1). 1000 did for every model except the two Ebbinghaus LNR fits, whose chains were still drifting through the sampling phase (4.10 C). The 2026-10-13 refit of everything runs 1500, the most the slowest models afford inside `long`'s 8 days: a warmup iteration costs 0.75-1.5x a sampling one, and 2000 would put the slowest LNR-6 / RDM-5 chains at 6.8-7.9 days (4.11) |
| draws | `IGC_SAMPLES=1000` (500 until 2026-10-09) | 8 chains x 1,000 = 8,000 draws per model. At 500, `b_ndt_Intercept` had a bulk ESS of only 53-310, and the per-participant solutions never switched; the longer run is the test of whether they ever do (4.10). It costs about a third more wall time: a `gam_lnr` shard goes from 43 h to about 57 h, and DDM-5 from a projected 2.5-5.2 days to 3.3-6.9, still inside `long`'s 8 days. It doubles the combine's memory, hence `combine.slurm`'s 320G (4.8). A task reads the default when it starts, so do not push a change to it while an array has pending tasks; `combine_model.R` refuses shards of different lengths |
| array | `--array=1-4` per model | 4 x 16 = 64 CPUs per model against `long`'s 140-CPU per-user cap, so 8 of the 12 production tasks run at once and the rest start as they finish (`--time` is per task, so pending costs nothing) |
| chains | `IGC_CHAINS=2`, 8 threads each | threading scales nearly linearly, so more chains per task buys nothing; one chain per task risks the whole task on one init (4.5) |
| memory | `--mem=32G` | ~6.5 GB per chain at full data (4.3), so ~13 GB used; the rest is headroom against an OOM killing a 20 h job |

`gam_ddm5` was not part of the original final run — it carries one more 2-D
smooth than DDM-4, and at 2.8x DDM-4 per gradient (4.7.1) it projects to
~2.5-5.2 days per shard against `gam_lnr`'s measured 42.7-44 h, which is close
enough to `long`'s 8-day wall to have wanted its own `--time`. **Submitted
2026-09-22 on `sussexneuro` instead** (job 11404322, `--array=1-4`, production
defaults otherwise), where the wall is 60 days and the question does not arise
— and where it runs in parallel with the two `long` arrays rather than behind
them (4.2.1). All four shards started immediately; account total 192 CPUs,
which the 140-CPU `long` cap could not have allowed.

Two things to watch on it, beyond the usual 6.1 checklist: the per-gradient
cost ratio at full data (4.7.1's 2.8x is a standalone density benchmark; the
30-participant end-to-end smoke said only ~1.4x, and open question 3 wants the
real figure), and the fact that its 64 CPUs come out of a **group** pool, so an
overrun inconveniences colleagues rather than only us. If a future DDM-5-sized
job needs bounding for that reason, `--time=8-00:00:00` is the exception to 3.7
worth making — the projection has 1.5-3x headroom inside 8 days.

`gam_ddm7`, `gam_rdm`, `gam_rdm5` and `gam_lba` (added 2026-09-18) are for the
**second account** — the CPU quota is per user, so a colleague fits those three while
`dmm56` fits the LNR pair and DDM-4. See README -> "Running from a second
cluster account". Since 2026-09-22 one of them can also come back to `dmm56` on
`sussexneuro`, whose quota is separate from `long`'s (4.2.1):
`./hpc fit <model> --partition=sussexneuro` runs as a further job alongside the
production arrays. Agree who takes which model either way — the point of
splitting was never the accounts, it was the quotas.

Two things about them that are decisions rather than defaults, and that would
be expensive to discover from the output:

- `gam_lba` uses `cogmod_lba2`. `cogmod_lba1` has one drift and no second
  accumulator, so it cannot model the choice in `dec(Error)`.
- `sigmazero = 1` in `gam_lba` is the LBA's **scaling constraint**, not a
  simplification: the likelihood is invariant to multiplying every drift, the
  start-point range and the boundary by a common constant, so exactly one
  parameter has to be pinned. `sigmaone` stays free. Freeing `sigmazero` too
  puts the chains on a ridge, which will read as poor Rhat rather than as an
  error.

`gam_ddm7` frees all three between-trial variabilities and smooths each. They
are identified through the shape of the RT distribution rather than its
location, weakly so even with flat predictors, and here each carries 25 tensor
coefficients plus a participant intercept. Read Rhat and ESS on `sigmadrift`,
`sigmabias` and `sigmandt` before trusting its `waic`; of everything in the
registry this is the one most likely to need reparameterising rather than more
hours.

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
  is incomplete at this scale — reopen the init investigation (4.5).
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
- **Not** many short chains on `general` to exploit its larger CPU allowance —
  warmup is per chain and would multiply. See 6.3, which works the arithmetic.

### 6.3 Many short chains on `general` instead of few long ones — no, not at full data

A recurring and reasonable-sounding proposal: `general` allows 400 CPUs against
`long`'s 140, so if a chain could manage a few hundred iterations inside the 8 h
wall, why not run far more chains of far fewer draws and let `combine_model.R`
pool them? The shards are combined anyway.

It does not work at full data, and the reason is worth writing down because the
idea will come up again.

**Warmup is per chain and cannot be amortised.** For `N` total draws at `d`
draws per chain the cost is

```
N x D  +  (N / d) x W
```

— sampling scales with the number of draws, warmup with the number of *chains*.
Halving `d` doubles the warmup bill. Going from 500 draws per chain to 125
costs **4x the warmup** to buy **2.86x** the CPUs (400 / 140). That is a net
loss before the wall is even considered.

**And the wall is binding anyway.** Measured on the first production run
(job 11385713, 2026-09-18): at full data `gam_lnr` had not printed iteration
100 after 1 h 15 m, i.e. **>40 s per warmup iteration**, matching the ~58 s of
the 4.4.1 pilot. Warmup 1000 alone is therefore **12-17 h** before a single
retained draw. No arrangement of chain counts fits that into 8 h, because every
chain pays it independently. The only mechanism that would amortise adaptation
across chains is carrying a metric and step size over from a pilot — a warm
start, closed by decision (5.5).

**The statistical objection is the more serious one.** 4.4.1 showed `lp__`
still rising monotonically after 300 warmup iterations at full data: those
chains had not reached the typical set. Pooling many such chains does not
average the transient away — it returns 4,000 draws of the same bias, and the
extra draws make the posterior look *tighter* while being wrong. Split-Rhat off
60-draw halves is too noisy to catch it. Whenever warmup dominates the chain
(here ~70% of it), few well-adapted chains beat many short ones.

**Where the idea is right:** below full data. At <= 480 participants a chain
finishes in about 1 h 48 m (4.4), so `general` plus many tasks is the better buy
— 6.6 already says this. It is also the right home for the `gam_ddm7`
reparameterisation experiments, where short exploratory chains are exactly what
is wanted.

**Free throughput that *is* available:** the per-partition quotas are separate.
On 2026-09-18 three `short` jobs started immediately while `gam_ddm4` sat in
`QOSMaxCpuPerUserLimit` against the `long` cap. So work that genuinely fits in
2 h or 8 h can run alongside the full-data jobs at no cost to them. The
constraint is the wall, not the quota.

Since 2026-09-22 there is a **better** version of that, and it has no wall to
work around: `sussexneuro` carries its own 256-CPU group quota and a 60-day
limit, so a full-data model can run there in parallel with the `long` arrays
rather than in the gaps around them (4.2.1). It does not change any of the
arithmetic above — warmup is still per chain, and many short chains are still
the wrong trade — it just means a *fourth* long-running job is free.

**The lever that does convert CPUs into a shorter wall** is threads per chain,
not more chains: 4.4.1 found within-chain threading scales nearly linearly,
while 4 chains x 4 threads bought no throughput over 2 x 8. That shortens a
chain instead of multiplying warmup. It still does not reach 8 h at full data
(halving ~19 h leaves ~10 h), and 1 chain per task stakes the task on one
chain's survival — defensible now that 0.3.3 has fixed the init failures (4.5),
but unproven at scale.

### 6.4 Ebbinghaus queue (2026-10-04)

cogmod 0.3.4 (`e01f8de`) was installed on 2026-10-04 with nothing pending on
the account, and the Ebbinghaus models still unfitted went out on it. The
`gam_lba*` models are left to the colleague who is fitting them. Two or three
models at a time is what `long` takes in practice (one `long` array is 64 CPUs
against the 140-CPU cap, and `long` is shared with other groups), so the
queue is limited to that, with the rest to be submitted as slots free up:

| model | state |
| --- | --- |
| `gam_ddm5_ebbinghaus` | done: 11423966 (task 4) + refit 11430275 (tasks 1-3, failed saves), combined and predicted 2026-10-09 |
| `gam_ddm4_ebbinghaus` | done: 11423962, combined and predicted 2026-10-07 (second combine, 11462276; the first failed at save) |
| `gam_rdm_ebbinghaus` | done: 11423963 + refit 11430274 (tasks 1 and 4, failed saves), combined and predicted 2026-10-09 |
| `gam_rdm5_ebbinghaus` | fitted: 11553382 (all 4 shards done 2026-10-12, 2-3.3 days each); combined (11595462, `--mem=192G`, 500 draws a chain, 30 min) and predicted (11595463) 2026-10-13. **Converged** (4.10) |
| `gam_lnr6_ebbinghaus` | 11553384 on `long` (zen5): all 4 shards done 2026-10-13 (3.4-3.9 days). **Not converged**: the `gam_lnr_ebbinghaus` adaptation failure (4.10 C), step sizes 1e-4 to 9e-6, Rhat 3.05, one chain frozen (shard 4 chain 2, lp__ SD 3.9). Combined anyway (11596092, 320G, 18 min), predicted (11596093) and priors (11596094) 2026-10-13, and in `2_modelcomparison.qmd` flagged NOT CONVERGED like the ML and VH LNR-6. Superseded by the 4.11 refit |
| `gam_lnr_ebbinghaus` | refit on 0.3.4: 11553388 on `sussexneuro` (zen5), combined 2026-10-12 (11594085), still not converged (Rhat ~3). Predicted anyway on 2026-10-13 (11595464), but left out of `2_modelcomparison.qmd`. The 0.3.3 shards, combined fit and predictions are under `..._20261002` |

Submit the two to-do ones as an array finishes, with
`./hpc fit <model> --constraint=amd_zen5`. Do not reinstall cogmod while any
array has pending tasks (2).

**Superseded 2026-10-13:** every model on every illusion is refitted on
cogmod 0.3.5, the 2026-10-12 data and 1500 + 1000 (4.11). The fits above are
the pre-refit generation; `combine_model.R` keeps the two apart by the
`m$cogmod` stamp and the draw count. They are now in `models/old_20261013/`
(3.6).

### 6.5 The 2026-10-13 refit: who runs what, and in what order

All 21 fits (the registry minus `gam_ddm7_mullerlyer`) are refitted on cogmod 0.3.5
(`eb4c6a3`), the 2026-10-12 data and 1500 warmup + 1000 draws (4.11).
Preparation on `dmm56`, 2026-10-13: the old generation moved to
`models/old_20261013/` (3.6), cogmod installed (2), and a 30-participant
smoke test of each of the seven families on `short`
(`smoke_refit_20261013/`). All seven compiled, started without an init
rejection, and saved: 4-17 min each, divergences 0 except LNR at 2% (warmup
300). The laptop's fits were deleted the same day. The four that existed only
there are in `old_20261013/laptop_only/`.

**The order follows what each fit can tell us, and how soon.**

- **`dmm56`'s `long`: the fast family, fastest first** (LNR, DDM-4, DDM-5, RDM:
  1.2-2.5 days per shard, 4.11). A problem with the new harness then shows up
  in 1-2 days rather than a week. `gam_lnr_verticalhorizontal`, the fastest
  fit in the registry, also tests the chain end to end: fit, combine, predict
  and priors on 0.3.5 and the new data.
- **`gam_lnr_ebbinghaus` in the first wave**, although it is the slowest of the
  fast family. It is the one failure that a longer warmup targets (4.10 C,
  open question 5), and `./hpc progress` shows its step sizes once warmup ends,
  after ~1.5 days.
- **`gam_rdm5_mullerlyer` at once, on `sussexneuro`.** It is the cheapest of
  the `sigmabias` models and the cleanest test of 4.10 B's open point: its
  MullerLyer fit split on 0.3.3, and the Ebbinghaus one converged on 0.3.4.
  If the full formula splits again, the simplified variants come next, and
  they are another week. So this test should not wait behind the fast family.
- **The other eight slow fits go to the two colleagues**: LBA x 3 plus
  RDM-5 Ebbinghaus to one, LNR-6 x 3 plus RDM-5 VerticalHorizontal to the
  other. Each submits all four at once. Their `long` cap runs two models at a
  time, so the second pair starts as the first finishes, ~10-12 days in all.
  The fourth model costs them no wall time: without it, half their slots
  would idle in the second wave. They need about four check-ins.
  Instructions: `refit_colleagues.md`.
- **zen5 is the binding resource, not the quotas.** Every long fit is pinned
  to the three rtx nodes (384 CPUs, shared with other users), because the
  slow family only fits inside 8 days there. With `dmm56` at 192 CPUs and
  two colleagues at 128 each, the demand exceeds the pool. Expect
  `PENDING (Resources)` on whichever account submits last.

| model | where | job | state |
| --- | --- | --- | --- |
| `gam_lnr_verticalhorizontal` | `dmm56` `long`, zen5 | 11596772 | running since 2026-10-13 18:58 (rtx-00) |
| `gam_lnr_ebbinghaus` | `dmm56` `long`, zen5 | 11596773 | running since 2026-10-13 18:58 (rtx-00, rtx-02) |
| `gam_rdm5_mullerlyer` | `dmm56` `sussexneuro`, zen5 | 11596823 | running since 2026-10-13 20:20 (rtx-02). First submitted as `gam_rdm5` (11596782, 18:59) and cancelled after 1.4 h for the rename |
| `gam_ddm4_ebbinghaus` | `dmm56` `long`, zen5 | 11596774 | pending behind the first two (`QOSMaxCpuPerUserLimit`) |
| `gam_ddm4_verticalhorizontal` | `dmm56` `long`, zen5 | 11596775 | pending |
| `gam_ddm5_verticalhorizontal`, `gam_lnr_mullerlyer`, `gam_rdm_ebbinghaus`, `gam_ddm5_ebbinghaus`, `gam_ddm5_mullerlyer`, `gam_ddm4_mullerlyer`, `gam_rdm_verticalhorizontal`, `gam_rdm_mullerlyer` | `dmm56` `long` (or `sussexneuro` once `gam_rdm5_mullerlyer` is done) | — | to submit in this order once the first `REPORT` lines are clean |
| `gam_lba_mullerlyer`, `gam_lba_ebbinghaus`, `gam_lba_verticalhorizontal`, `gam_rdm5_ebbinghaus` | colleague A | — | instructions written |
| `gam_lnr6_mullerlyer`, `gam_lnr6_ebbinghaus`, `gam_rdm5_verticalhorizontal`, `gam_lnr6_verticalhorizontal` | colleague B | — | instructions written |

Update the table as things move, and add a `node_history.md` row for every
shard that finishes. The colleagues send each model's combined fit,
prediction file, priors file and `./hpc log` output. Upload the fit to
`dmm56`'s `models/combined/` (md5-checked) so that `./hpc priors` and
`mixing_participants.R` see every model in one place.

---

## 7. Open questions

1. **How long is the full-data cold-start transient, and does warmup 1000
   clear it?** The warmup-300 pilot (4.4.1) shows > 300 iterations at 1023
   leapfrog steps (~80 s each) and no adaptation. Warmup 1000 has never been
   run at full data. The first production `gam_lnr` job answers this — see 6.1.
   *Interim, 2026-09-18:* 1 h 15 m in, job 11385713 had not reached iteration
   100, i.e. >40 s per warmup iteration, so warmup 1000 costs 12-17 h. That is
   consistent with the 12-20 h per chain projected in 4.4.1 and comfortably
   inside the 2-day `--time`, but it is also what rules out the 8 h partition
   for full data (6.3). *Answered (4.11):* warmup 1000 adapted every model
   except the Ebbinghaus LNR pair, and costs 10-56 h per chain depending on
   the model.
5. **Does warmup 1500 clear the Ebbinghaus LNR drift?** (4.10 C, 4.11) Both
   LNR models on that illusion ended warmup 1000 with step sizes ~1e-4 and
   chains still moving. If 1500 does not fix it, reparametrise rather than
   lengthen again.
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
| `hpc` | the driver — check/setup/push/install/precompile/models/fit/combine/predict/priors/archive/progress/queue/log/ls/pull/cancel/sh |
| `refit_colleagues.md` | step-by-step for a colleague's account in the 2026-10-13 refit (6.5) |
| `install_pkgs.R` | builds the project R library (`./hpc install`); enforces cogmod >= 0.3.5 (`IGC_COGMOD_MIN`) |
| `precompile.R` | builds the CmdStan precompiled header (`./hpc precompile`) |
| `models.R` | **the model registry** — one entry per model, read by everything else |
| `fit_model.R` | fits the model named by `IGC_MODEL`, one shard per array task |
| `fit.slurm` | array job for the above; its `#SBATCH` defaults are the production config |
| `combine_model.R` | merges one model's shards, adds `loo`, keeps them |
| `combine.slurm` | job for the above |
| `README.md` | command reference, partition quotas, PCH detail |
| `keep_stored_basis.R` | `keep_stored_basis()`, the local workaround for 3.9 |
| `hpc.local` | **gitignored** — this machine's account settings |
| `server.md` | **gitignored** — local notes on this project's cluster directories (account, keys and VPN live in the lab hub's `hpc/private/`) |
