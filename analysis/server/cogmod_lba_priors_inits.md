# cogmod: LBA drift priors and starting values to investigate and fix

2026-10-01. This is a to-do list for cogmod, prompted by the full-data `gam_lba`
fit splitting into two posterior modes (`cogmod_lba_modes_issue.md`). That
investigation found that **cogmod's LBA2 Stan code is accurate** (to 1e-11)
and is not the cause. The cause is the model's geometry. The defaults below did
not create the modes either, but they are where cogmod either works against
the data or makes an arbitrary choice, and they are cheap to change.

**The plan:** adjust these in cogmod, then re-run `gam_lba` to see whether the
split replicates. If it does, the fallback is a model change, not a cogmod one
(see the end).

**Status 2026-10-02: cogmod 0.3.4 (`dev`, `e8c0555`) has made these
changes.** See its `NEWS.md`.

| item | in 0.3.4 |
| --- | --- |
| A1 `driftone` start | **done**, option 1: `driftone = 1` |
| A2 `sds` start | **done**, option 1: 0.1 link units of wiggle |
| A3 jitter | **changed after all**: the jitter on regression blocks is divided by their design's largest row norm. On these models the `bs_*` start had tilted every smoothed dpar by a median of 2.3 link units. |
| B1 `sds` prior | **done**, option 1: `exponential(rate)` per smooth, set from its basis; it now applies to every smooth, `mu`'s included |
| B2 `mu` / `driftone` asymmetry | **done**, option 1: `mu`'s participant SD and `sds` match `driftone`'s; intercept and slopes are still left to brms |
| B3 | **docs done** (`?rcogmod_lba2` opens with the `sigmaone ~ 1` advice); the warning is not implemented |
| (new) unpenalised smooth columns | `bs_*` on a family-named dpar get the dpar's slope prior in link units |

What is left is §D.3, the real re-run. 0.3.4 changes the priors of every model
with a smooth, not only the LBA's, so read `AGENT.md` §2 before
`./hpc install cogmod`.

The scripts behind the numbers are in `lba_modes/` and run from
`analysis/server/`. cogmod references are to the `dev` branch at `fc0e9f5`;
the cluster runs `d04c7f8`, and the code below is the same in both.

---

## A. Starting values

### A1. LBA2 starts `driftone` equal to `mu`: a 50/50 race

**Where:** `R/core_choice.R`, the `cogmod_lba2` registry entry:
`init = list(mu = 3, driftone = 3, sigmazero = 1, sigmaone = 1, sigmabias = 0.5, boundary = 0.5)`.

**Why it is suspect.** The RDM entry in the same file was changed to
`driftone = 1`, "a third of `mu`", after the RDM brittleness work.
`?cogmod_inits` gives the reason: starting the error accumulator too fast
"costs hundreds" of log-density units, which a cold chain turns into momentum
along the flat `driftone` direction until its step size collapses. Too slow
costs a few dozen. The LBA2 never got the same change.

The error rate implied by the LBA2 start, with the other dpars at their inits
(`init_error_rate.R`):

| `driftone` start | 3 (now) | 2 | 1.5 | 1 | 0.5 | 0 | -1 | -2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| implied error rate | **50%** | 27% | 18% | 11% | 6.9% | 4.2% | 1.7% | 0.8% |

Few real 2AFC data sets sit at 50% errors. This one is at 23.7% overall,
ranging from 1.2% to 81% across cells.

**Options:**
1. **Minimal:** `driftone = 1`, matching the RDM (11% errors at start).
2. **Data-informed:** solve for the `driftone` that reproduces the observed
   share of `dec == 1` at the other start values, as `.ndt_start()` already
   takes `ndt` from the data. `cogmod_inits()` has `make_standata()` in hand,
   so `dec` is available. Values: 1% → -1.66, 5% → 0.18, 10% → 0.88,
   23.7% → **1.84**, 40% → 2.58. Clamp to a sane range, and lower `mu` instead
   when `dec == 0` is the rarer response.

**Caveat.** The note dated 2026-09-22 in `R/cogmod_inits.R` (above
`.ndt_start()`) records that moment-matching *every* start to the data bought
nothing. A start merely far from the mode is absorbed by the first metric
window; only a start *on a flat region* is not. 50% errors is not on the
plateau, so the gain may be modest. The case for the change is consistency
with the RDM and the asymmetry above (starting too fast costs more than too
slow). In `gam_lba` all 8 chains started within jitter of the same point and
still split, so this alone is not expected to fix it. It needs testing (§D).

### A2. Smooth SDs start at 0.05, which is "flat" in a way that depends on the basis

**Where:** `R/cogmod_inits.R`, in `.init_value()`:
`if (startsWith(d$name, "sds_")) return(rep(0.05, n))`. The comment
justifies 0.05 by "tensor-product basis values in the tens".

**What the basis actually looks like here.** For this `t2(..., k = c(5, 5),
bs = c("cr", "cr"))` smooth, one unit of `sds` moves the linear predictor by
**0.08–0.12 link units** (RMS over rows of `sqrt(rowSums(Zs^2))`,
`smooth_scale.R`). That barely depends on data size (`smooth_scale_n.R`):

| participants | 30 | 120 | 480 | 2,215 |
| --- | --- | --- | --- | --- |
| rows | 3,810 | 15,100 | 60,600 | 324,000 |
| link units per unit `sds` (penalties 1 / 2 / 3) | 0.118 / 0.097 / 0.097 | 0.118 / 0.097 / 0.097 | 0.107 / 0.087 / 0.087 | 0.096 / 0.078 / 0.079 |

So a start of 0.05 is ~0.005 link units of wiggle: flat. The posterior `sds`
of the drift smooths' largest penalty is 10–37 (table in B1), 5.3–6.6 log
units above the start. Starting flat is a defensible choice. But "values in the tens" does not
describe this basis, so the constant was calibrated on a different one.

**Options:**
1. Express the start in link units: read the `Zs_*` scale from the
   `make_standata()` output `cogmod_inits()` already builds, and set
   `sds = target / scale`, e.g. a target of 0.1 link units.
2. Find which model the "tens" measurement came from (an `s()` with large `k`?
   a `te()`?), and check whether a single scale-aware rule covers both.

### A3. Jitter: no change

The population tier (0.25) is what makes Rhat meaningful, and it was not the
issue. Chains starting within it still split.

---

## B. Priors on the drift smooths

### B1. `exponential(1)` on `sds` is a tight prior in link units for this basis

**Where:** `R/cogmod_priors.R`, `.priors_dpars()`. The 0.3.3 fix replaces
brms's blanket `sds` row with the family's `exponential(1)` for every
non-response dpar.

**What it means here.** With the basis scale of A2, the prior medians
translate to these amounts of wiggle per t2 penalty:

| prior on `sds` | median | link-unit wiggle (SD over rows) |
| --- | --- | --- |
| `exponential(1)` (cogmod, every dpar but `mu`) | 0.69 | **0.05–0.07** |
| `student_t(3, 0, 2.5)` (brms, kept on `mu`) | 1.9 | 0.15–0.18 |

**What the posterior used** (`smooth_scale.R`, the largest t2 penalty, posterior
means, modes A / B):

| dpar | link | posterior `sds` | wiggle (link units) | whole-smooth range |
| --- | --- | --- | --- | --- |
| `mu` | identity | 13.4 / 10.5 | 0.64 / 0.56 | 4.6 / 4.3 |
| `driftone` | identity | 14.7 / 36.6 | **1.50 / 5.10** | **23.7 / 27.5** |
| `sigmaone` | softplus | 4.4 / 8.2 | 0.26 / 0.50 | 2.0 / 1.8 |
| `sigmabias` | softplus | 4.1 / 13.0 | 0.30 / 1.16 | 2.3 / 6.3 |
| `boundary` | softplus | 4.7 / 5.0 | 0.16 / 0.34 | 2.4 / 1.8 |
| `ndt` | log | 2.4 / 1.6 | 0.11 / 0.09 | 1.9 / 0.9 |

Against `exponential(1)`'s 0.05–0.07, `driftone`'s smooth uses 25–85x the
prior's scale and `sigmabias`'s 5–19x. Only `ndt`, on a log link, is close
(1.5–2x). `mu`, under brms's wider half-t, uses 3–4x. At 324k rows the likelihood wins regardless, but mode B
pays ~50 units of prior for it, and the prior is plainly not "weakly
informative" on the drift scale. At a few thousand rows it would bind.

The `exponential(1)` choice was argued for log- and logit-linked dpars (a
half-t median of 1.9 "lets the smooth alone walk a `sigmabias` across its
whole range"). That argument does not carry over to identity-linked drifts,
and it was made without the basis scale in view.

**Options:**
1. **Scale-aware rate:** set the `sds` prior per smooth term in link units, by
   dividing by the basis scale read from `make_standata()`. brms accepts a
   prior per `sds` row via `coef`.
2. **Link-aware default:** keep `exponential(1)` for log/logit/softplus dpars
   and widen it for identity-linked drifts (`driftone`; `mu` if B2 is
   adopted).
3. Either way, check with a prior predictive: draws of the smooth on the data's
   covariates should span a plausible range of drift (a few units), not 0.05
   or 50.

### B2. The two drifts get different priors on every class

`mu` (the correct accumulator's drift) and `driftone` (the error
accumulator's) are the same kind of quantity on the same scale. In this fit
(`priors_used.R`):

| class | `mu` | `driftone` |
| --- | --- | --- |
| Intercept | `student_t(3, 0.6, 2.5)` (brms) | `normal(1, 2)` (cogmod) |
| `b`, incl. the smooth's linear part | flat in cogmod; `normal(0, 1)` added by our `fit_model.R` | `normal(0, 1.5)` (cogmod) |
| `sd` (Participant) | `student_t(3, 0, 2.5)` (brms) | `exponential(1)` (cogmod) |
| `sds` | `student_t(3, 0, 2.5)` (brms) | `exponential(1)` (cogmod) |

This follows a deliberate rule in `.priors_dpars()`: the response's own
predictor (`dpar == ""`) is left to brms. For a race family, though, `mu` is
not "the response". It is accumulator 0's drift, symmetric with accumulator
1's, and which accumulator is `mu` is just the `dec()` coding. Swap the coding
and every prior on the two drifts swaps with it.

**Options:**
1. For the race families (`cogmod_lnr`, `cogmod_rdm`, `cogmod_lba2`), treat
   `mu`'s `sd` and `sds` like the other accumulator's. Leave slopes alone if
   the "don't shrink the effects of interest" principle is to be kept.
2. At minimum, document the asymmetry in `?cogmod_priors`. The docs already
   say to mirror `driftone`'s prior onto `mu` by hand when the rare response
   is on `mu`.

### B3. `driftone`'s intercept and slopes sit 3–4 prior SDs out

| | prior | posterior (A / B) |
| --- | --- | --- |
| `driftone` intercept | `normal(1, 2)` | -4.9 / -7.0 |
| first linear smooth coefficient (`bs`, Xs column SD 0.8) | `normal(0, 1.5)` | 5.0 / 4.7–5.0 |
| `sigmaone` intercept (softplus) | `normal(0, 1)` | 2.0 / 2.6 |

These are not a prior bug. They are where the truncated-normal plateau puts an
accumulator that rarely wins once its SD is free (`?rcogmod_lba2`, "Negative
drift rates"). Here that happens on the facilitating trials (1–3% errors),
where the error accumulator sits at v/s of -4 to -5.4. Any design with an easy
condition will do this when `sigmaone` has predictors. A tighter prior would
only bias the estimate. What cogmod could add instead:

- **A warning** along the lines of `.warn_scale_ray()`: when the free drift SD
  gets predictors, warn that wherever that accumulator rarely wins, its drift
  and SD are identified only through |v|/s², and suggest `sigmaone ~ 1` or a
  fixed `sigmaone`.
- **Docs:** make that the headline advice in `?rcogmod_lba2`, not a closing
  remark.

---

## C. Not to change

- **The Stan density.** It matches direct integration to 9.7e-12 across both
  modes' parameters and RTs of 0.05–3.5 s, and its gradients are smooth through
  the drift's sign change (`lba_numerics.R`). rtdists is the less accurate of
  the two at very short RTs.
- **The truncation of drifts at zero.** It is the field's convention, and
  without it the density does not integrate to one.
- **The `b = boundary + sigmabias` parametrisation.** The ridge it leaves is
  documented and inherent to the LBA.

---

## D. How to test, cheapest first

1. **Prior predictive** for B1: draw smooths from the old and new `sds`
   priors on the `t2` basis and report the spread on the drift scale. Minutes.
2. **Simulation benchmark** for A1, A2, B1 and B2, in the style of
   `benchmarks/inits_ablation/`. Simulate LBA2 data whose error rate spans
   ~1% to ~80% across a 2-D design like this one, about 100–200 participants
   x ~150 trials. Fit with the old and the new defaults, 8 chains each, and
   record:
   - whether chains split (spread of lp__ between chains, population Rhat);
   - warmup leapfrogs and adapted step size;
   - parameter recovery, especially `driftone` / `sigmaone` where errors are
     rare.

   This is the real test of A1 and A2. The caveat in A1 predicts no change,
   and that prediction is worth checking.
3. **The real re-run:**
   - `./hpc install cogmod` to pull the change onto the cluster.
   - Run `gam_lba` with its own `IGC_MODELS_DIR` or `IGC_FILE_REFIT=always`.
     `file_refit = "never"` silently reuses a shard already on disk.
   - Compare chain agreement against the 2026-10-01 fit: the population Rhat
     (13.8 then), lp__ per chain, the `driftone` intercept.

**If the split replicates:** the next step is the model change kept in reserve,
`sigmaone ~ 1 + (1 | Participant)` (`cogmod_lba_modes_issue.md` §6, option
a). With a drift SD shared across conditions, |v|/s² pins `driftone` in every
condition, which removes the condition-by-condition flat direction the modes
use.

---

## Appendix: the model as fitted

**Data.** MullerLyer only: 323,981 trials, 2,215 participants. `RT` is in
seconds. `Error` is 0/1; `dec(Error)` makes 0 (correct) accumulator 0, whose
drift is `mu`, and 1 (error) accumulator 1, `driftone`. Both predictors are
normalised within illusion (`fit_model.R`):
`Illusion_DifferenceZ = 2 * normalize(|Illusion_Difference|) - 1`, from -1
(hardest) to 1 (easiest), and
`Illusion_StrengthZ = sign(Illusion_Strength) * normalize(|Illusion_Strength|)`,
from -1 to 1, with 0 for no illusion and > 0 for conflicting.

**Formula** (`models.R`, `gam_lba`, written out):

```r
f <- brms::bf(
  RT | dec(Error) ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
  driftone        ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
  sigmaone        ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
  sigmabias       ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
  boundary        ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
  ndt             ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
  sigmazero = 1,
  poutlier        ~ 1 + (1 | Participant),
  family = cogmod::cogmod_lba2()
)
```

The links are cogmod's defaults: `mu` and `driftone` identity; `sigmaone`,
`sigmabias` and `boundary` softplus; `ndt` log; `poutlier` logit. The
`(1 | Participant)` terms are separate, so the participant intercepts are
uncorrelated across dpars.

**The process.** On each trial, accumulator 0 (correct) draws a drift
v0 ~ N+(`mu`, 1) and accumulator 1 (error) draws v1 ~ N+(`driftone`,
`sigmaone`), each a Normal truncated at zero. Each starts at its own
z ~ U(0, `sigmabias`) and rises linearly to b = `boundary` + `sigmabias`. The
response is the first to arrive, at RT = `ndt` + min_k (b - z_k) / v_k. With
probability `poutlier` the trial is instead a guess: a uniform choice, with a
half-Normal RT of scale 0.2 s.

**Priors.** `cogmod_priors(f, data)`, plus `normal(0, 1)` on `mu`'s slopes,
which cogmod leaves flat (`fit_model.R`). As resolved per class:

| dpar | Intercept | `b` (incl. smooth linear part) | `sd(Participant)` | `sds` |
| --- | --- | --- | --- | --- |
| `mu` | `student_t(3, 0.6, 2.5)` | `normal(0, 1)` | `student_t(3, 0, 2.5)` | `student_t(3, 0, 2.5)` |
| `driftone` | `normal(1, 2)` | `normal(0, 1.5)` | `exponential(1)` | `exponential(1)` |
| `sigmaone` | `normal(0, 1)` | `normal(0, 0.5)` | `exponential(1)` | `exponential(1)` |
| `sigmabias` | `normal(0, 1)` | `normal(0, 0.5)` | `exponential(1)` | `exponential(1)` |
| `boundary` | `normal(0, 1)` | `normal(0, 0.5)` | `exponential(1)` | `exponential(1)` |
| `ndt` | `normal(-1.2, 0.5)` | `normal(0, 0.2)` | `exponential(1)` | `exponential(1)` |
| `poutlier` | `normal(-5, 1)` | — | `exponential(1)` | — |

**Inits and sampler.**
- Inits: `cogmod_inits(f, data)` with the default jitter (0.25 on
  population-level parameters, 0.05 on the hierarchical ones). The targets are
  `mu = 3, driftone = 3, sigmaone = 1, sigmabias = 0.5, boundary = 0.5,
  poutlier = 0.02`, with `ndt` at half the first RT percentile and `sds` at
  0.05.
- Sampler: `brm(..., stanvars = cogmod_stanvars(f), backend = "cmdstanr")`,
  warmup 1000 + 500 draws, 2 chains x 8 threads per shard, 4 shards. Stan
  defaults otherwise: `adapt_delta = 0.8`, `max_treedepth = 10`, diagonal
  metric. Compiled with `stanc` O1 and `STAN_NO_RANGE_CHECKS`.

| script (`lba_modes/`) | what |
| --- | --- |
| `init_error_rate.R` | error rate implied by the LBA2 start; data-informed `driftone` (A1) |
| `smooth_scale.R` | link units per unit `sds`, and the posterior smooths in those units (A2, B1) |
| `smooth_scale_n.R` | basis scale against data size (A2) |
| `priors_used.R` | every prior the fit resolved (B2, appendix) |
| `lba_numerics.R` | Stan density against integration, and gradients (C) |
