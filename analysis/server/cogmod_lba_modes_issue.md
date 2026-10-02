# `gam_lba` has two posterior modes: what causes them, and what cogmod could change

2026-10-01. The fit is the combined `gam_lba_MullerLyer.rds`: 4 shards x 2
chains, warmup 1000 + 500, cold starts, brms 2.21.0, cogmod 0.3.3 (`d04c7f8`).
It was fitted on the second account. Status: **not converged, not predicted,
not in the model comparison** (AGENT.md §4.9). The scripts behind every number
here are in `lba_modes/` and run from `analysis/server/` on a laptop. None of
them refits anything.

## Summary

- **The two modes are real.** Mode A (5 chains) and mode B (3 chains) are
  separated by a valley at least 4,500 log-likelihood units deep along the
  straight line between them. Mode A is a tight, well-mixed *local* optimum
  (max Rhat 1.03 across its 5 chains). B fits the data better by
  **1,364 log-likelihood units** and is broad: its own 3 chains disagree (Rhat
  1.86).
- **cogmod's Stan code is not at fault.** The LBA2 density matches direct
  numerical integration to 1e-11 at both modes' parameter values, and its
  gradients are smooth everywhere a chain goes.
- **The model's geometry is the problem.** On facilitating trials the error
  accumulator almost never wins (1–3% errors), so it sits deep in the
  truncated-normal tail (v/s of -4 to -5.4). There its drift and drift SD are
  identified only through |v|/s². The two modes put that accumulator at
  different points of the ray (drift -7 to -12 against -15 to -17) with the
  same realised behaviour. They also split threshold and start-point range
  differently across the design.
- **Both modes predict almost the same behaviour.** In the conflicting and
  facilitating cells, their error rates agree to within 0.01, their medians to
  within 0.035 s and their 90th percentiles to within ~0.1 s. Both miss the
  same thing, by more than that. Longer runs would not reveal the "correct" parameters. The parameters
  that differ are the weakly identified ones.
- **The priors and inits are not the cause, but two are worth fixing.** At
  324k rows the data overwhelm every prior. The LBA2's default init for the
  error accumulator (`driftone = 3`, equal to `mu`) is the cold-start pattern
  cogmod already corrected for the RDM, and the LBA2 never got that fix.

## 1. The two modes are real

**The likelihood, not the prior, separates them.** `lba_decompose.R` splits
each chain's lp__ into the log-likelihood, `lprior`, the std-normal terms on
the non-centred group effects and smooth coefficients, and the Jacobians. z and
zs are rebuilt as r / sd and s / sds. Means over the 500 draws of each chain:

| mode | chains | lp__ | lprior | log-likelihood (up to a constant) |
| --- | --- | --- | --- | --- |
| A | 1, 2, 4, 6, 8 | 36,965 | -91 | 44,882 |
| B | 3, 5, 7 | 38,260 | -143 | 46,244 |

The A chains are within 20 units of one another. The B chains span 53 units
(46,221–46,274).

**A deep valley lies between them.** `lba_path.R` interpolates linearly, on the
unconstrained scale the sampler moves on, from the last draw of an A chain to
the last draw of a B chain. Both pairs tried show a deep dip:

| w (0 = A, 1 = B) | 0 | 0.3 | 0.5 | 0.7 | 1 |
| --- | --- | --- | --- | --- | --- |
| log-lik, chain 1 → 5 | 59,313 | 55,971 | 54,841 | 57,309 | 60,606 |
| log-lik, chain 6 → 3 | 59,316 | 55,016 | 53,619 | 56,500 | 60,735 |

A straight line overstates the barrier, but one this deep means chains in A are
not just slow along a ridge.

**A mixes and B does not.** Max Rhat over the population-level parameters is
1.03 across the 5 A chains and 1.86 across the 3 B chains. The B chains
disagree on the `ndt`, `sigmabias` and `boundary` intercepts and on the
`sigmabias` smooth. All 8 chains started within the population-level jitter
(SD 0.25) of the same point. They split during warmup, so with more chains both
basins would turn up again.

## 2. What mode B fits better

`lba_perobs.R` computes each observation's log-likelihood by hand from the
draws and the `standata()` design matrices, through `dcogmod_lba2()`. It does
not use brms's predict path, which the laptop's brms 2.23.1 gets wrong for this
fit. Validation: the per-draw sums match the lp__-implied log-likelihood for
all 8 chains, with the same constant offset (14,363.928) to the third decimal.

| | B − A (log-lik) |
| --- | --- |
| correct responses | **+1,433** |
| errors | -69 |
| conflicting D2–D4 | +990 |
| facilitating D1–D4 | +372 |
| conflicting D1 (80% errors) | -13 |

D1–D4 are quartiles of illusion difference, from hardest to easiest. Among
correct responses, B gains in 0.5–0.6 s (+1,397) and above 1.2 s (+1,475), and
loses below 0.4 s (-1,297). So B buys a better shape of the
**correct-response** RT distribution. The gain is not concentrated in a few
participants: the top 10% hold 32% of the |B − A| total, and 62% of
participants favour B.

## 3. Where the modes differ

Posterior means by cell (`lba_perobs_summary.R`):

| cell | errors | `driftone` A / B | `sigmaone` A / B | v/s A / B | realised mean drift A / B | `sigmabias` A / B | `boundary` A / B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| conflicting D1 | 81% | 2.8 / 3.5 | 1.7 / 2.2 | 1.8 / 1.7 | 3.2 / 4.1 | 1.06 / 1.76 | 0.83 / 0.71 |
| conflicting D3 | 30% | -0.2 / -1.8 | 1.8 / 2.3 | 0.0 / -0.5 | 1.37 / 1.40 | 0.89 / 1.37 | 1.14 / 0.92 |
| facilitating D1 | 2.7% | -6.9 / -15.5 | 1.9 / 2.9 | -3.7 / -5.4 | 0.48 / 0.53 | 0.61 / 0.41 | 0.77 / 0.99 |
| facilitating D4 | 1.2% | -12.0 / -17.3 | 2.6 / 3.4 | -4.7 / -5.0 | 0.51 / 0.63 | 0.60 / 0.13 | 0.73 / 1.30 |

- **The error accumulator is on cogmod's documented plateau.** Where it rarely
  wins, a Normal truncated at zero converges to an Exponential as its location
  runs to minus infinity with |v|/s² held fixed (`?rcogmod_lba2`, "Negative
  drift rates"). At v/s of -3.7 to -5.4, P(v > 0) is 1e-4 to 3e-8, so only the
  extreme tail is in play. The modes disagree by 5–9 units of drift there, yet
  the realised mean drift is nearly identical (0.48–0.63). That is the ray.
- **The smooths have to span both regimes.** The `driftone` surface runs from
  +3.5 to -17 across the design, because the "error" accumulator is the
  dominant response on hard conflicting trials and a near-silent one on
  facilitating trials. Its tensor smooth SD is 15 in A and 36 in B.
- **Start-point range and threshold are allocated differently.** B shrinks
  `sigmabias` toward 0.13–0.4 on facilitating trials and grows it to 1.0–1.8
  on conflicting ones, with `boundary` moving the other way. A keeps
  `sigmabias` at 0.6–1.06 everywhere. This is the threshold ridge cogmod
  documents (`b = boundary + sigmabias`), and it is where the correct-RT gain
  of §2 comes from.

## 4. Both modes predict the same behaviour

Simulated from the fit (4 draws per chain: 20 data sets for A, 12 for B).
Quantiles in the "none" cells, which have only 12–200 errors, are noisy and
left out:

| cell | | error rate | correct median | error 10% | error median | error 90% |
| --- | --- | --- | --- | --- | --- | --- |
| conflicting D2 | observed | 0.575 | 0.845 | 0.455 | 0.657 | 1.365 |
| | A / B | 0.581 / 0.572 | 0.821 / 0.824 | 0.445 / 0.441 | 0.672 / 0.671 | 1.374 / 1.380 |
| facilitating D3 | observed | 0.012 | 0.546 | 0.283 | 0.462 | 1.053 |
| | A / B | 0.013 / 0.014 | 0.546 / 0.547 | 0.121 / 0.120 | 0.609 / 0.581 | 1.419 / 1.384 |

Both modes misfit errors in the easy cells the same way. Their predicted errors
are far too dispersed (10th percentile 0.12 s against 0.28 s observed, 90th
1.4 s against 1.05 s). That is the Exponential-like drift distribution of an
accumulator on the plateau, which produces both very fast and very slow
finishes. Errors on conflicting D4 are too slow in both (median 0.69–0.71 s
against 0.58 s).

For context, the pooled 8-chain elpd_loo is 51,806, against 52,559 for
`gam_lnr6` (2-shard file) and 52,528 for `gam_ddm5`. The mode decides where the
LBA ranks.

## 5. The three questions

### 5.1 Is cogmod's LBA Stan code brittle?

Not numerically. `lba_numerics.R` compiles cogmod's own Stan `lpdf` and
compares it with direct numerical integration over the start point (and with
`rtdists::n1PDF()`). It covers both modes' parameters, both responses, and
decision times from 0.05 to 3.5 s. The largest difference is **9.7e-12**;
rtdists is the less accurate of the two, off by 0.1–0.3 at t ≤ 0.05 s.
Central-difference gradients along `driftone`, from +3 down to -20 (through
the sign change where `lsurv_trunc` switches branch), agree between steps of
1e-4 and 1e-6 to 5e-5 relative or better. The one exception is 7e-4, at a
point where the gradient itself is 5e-6. There are no kinks, no noise and no
cliffs.

The brittleness is in the parametrisation, which the code implements
correctly:

- **The truncation at zero** (the field's `posdrift = TRUE` convention)
  creates the ray for any accumulator that rarely wins. Dropping it would make
  the density integrate to less than one. That is a different model, not a
  fix, and it is not recommended.
- **The `b = boundary + sigmabias` parametrisation** leaves the threshold
  ridge. That is documented and inherent.

The only code-level change that would remove the ray is a different drift
distribution. cogmod already has one: `cogmod_lnr()` with a free `sigmabias`
is a ballistic accumulator with a Uniform start point and a **log-normal**
rate (core_shifted.R, "LogNormal accumulator with a start-point range"). That
is `gam_lnr6`, which is already fitted and leads the comparison.

### 5.2 Should the default priors and inits change?

**Priors: they are overwhelmed here, but they do not cause the modes.** The
prior difference between modes is 52 units, against 1,364 for the likelihood.
Where the posterior sits relative to cogmod's defaults:

| parameter | default prior | posterior (A / B) |
| --- | --- | --- |
| `driftone` intercept | `normal(1, 2)` | -4.9 / -7.0 (3–4 SD below) |
| `sigmaone` intercept (softplus) | `normal(0, 1)` | 2.0 / 2.6 |
| smooth SD, `driftone` (largest t2 term) | `exponential(1)` | 15 / 36 |
| smooth SD, `sigmabias` (largest t2 term) | `exponential(1)` | 4 / 12–14 |

For a typical application of a few thousand trials these defaults are sensible
and would regularise. At 324k rows no weakly informative default can steer
geometry like this. Only a parametrisation can. One default is worth a look in
cogmod: `exponential(1)` on `sds` for identity-link drift smooths. Here it
sits 15–36 prior means below the posterior. That is harmless for convergence,
but it is not weakly informative on that scale.

**Inits: one concrete fix, though it is unlikely to solve this fit on its
own.** The LBA2 registry entry (core_choice.R) starts `mu = 3, driftone = 3,
sigmazero = 1, sigmaone = 1`, a 50/50 race everywhere. cogmod's RDM entry
starts `driftone` at a third of `mu` for exactly this reason. `?cogmod_inits`
says too fast "costs hundreds" of log-density units, which a cold chain turns
into momentum along the flat `driftone` direction until its step size
collapses. The LBA2 never got that change. A data-informed start (the observed
error rate, as `ndt` is taken from the data) would be the general version. But
every chain here started at the same point and still split, and the posterior
`driftone` spans +3.5 to -17 across cells. A different start would most likely
change *which* basin chains land in, not remove the valley.

### 5.3 Would longer or subset runs help?

Not to learn the "correct" parameters. The parameters that differ between the
modes are the error accumulator's drift and SD on the plateau, and the split of
threshold against start-point range. They are weakly identified by
construction, both modes reproduce the behaviour equally (§4), and a subset
identifies them less, not more. More warmup at full data costs 2.3–3 days a
shard and would add chains to whichever basin they fall into.

A subset run is useful for something narrower: as a cheap **screen** of a
reparametrised variant. If its chains split at ~200 participants, the variant
is rejected. If they agree there, that is necessary but not sufficient for full
data.

## 6. Options

For cogmod, independent of this project:

1. **Inits:** start LBA2's `driftone` below `mu`, as the RDM does, or from the
   observed error rate. Cheap and general.
2. **A warning:** when the non-pinned drift SD (`sigmaone`, or `sigmazero` if
   `sigmaone` is the pinned one) gets predictors, warn that wherever that
   accumulator rarely wins, its drift and SD are identified only through
   |v|/s². This follows the pattern of the existing `.warn_scale_ray()`.
3. **The `sds` prior:** revisit `exponential(1)` for identity-link drift
   smooths. Flag rather than change. It is not what broke this fit.

For this project, the LBA variants that remove the freedom the modes use, from
least to most restrictive:

- **(a) `sigmaone ~ 1 + (1 | Participant)`.** A drift SD shared across cells
  pins `driftone` through |v|/s² in each cell, which removes the cell-wise
  ray. The accumulators can still differ in SD.
- **(b) Also `sigmabias ~ 1 + (1 | Participant)`.** This removes the
  threshold/start-point surfaces the B chains disagree on.
- **(c) `sigmaone = 1`.** The textbook single-`sv` LBA. It removes the ray
  entirely, at the cost of the error-accumulator SD the data clearly want
  (2–3).
- **(d) Drop `gam_lba` from the comparison.** Say why, and point to
  `gam_lnr6` as the log-normal-rate ballistic accumulator.

**Decided 2026-10-01:** cogmod's defaults are adjusted first
(`cogmod_lba_priors_inits.md`), and `gam_lba` is re-run unchanged to see
whether the split replicates. Option (a) is the next step only if it does.

Each variant is a new `models.R` entry. The screen would be about 200
participants with 4 chains x 8 threads on `long`; scaling the measured
full-data cost puts that at about 4–7 h per chain. A full-data fit after that
would be 4 shards x 2–3 days.

## Reproduction

From `analysis/server/`, with the combined fit in `analysis/models/`:

| script | what | time |
| --- | --- | --- |
| `lba_modes/lba_decompose.R` | lp__ split per chain (§1) | ~1 min |
| `lba_modes/lba_path.R` | the A→B path (§1) | ~3 min |
| `lba_modes/lba_perobs.R` | per-observation log-lik, cell parameters, simulations (§2–4) | ~4 min |
| `lba_modes/lba_perobs_summary.R` | the tables of §2–4 | ~1 min |
| `lba_modes/lba_coefs.R` | smooth coefficients and SDs per chain, within-mode Rhat | ~1 min |
| `lba_modes/lba_numerics.R` | Stan density against integration and rtdists, gradients (§5.1) | ~2 min, compiles Stan |
| `lba_modes/data_desc.R` | error rates and RTs by cell | seconds |

The output `.rds` files land in `lba_modes/` and are gitignored. The
per-observation one is 128 MB.
