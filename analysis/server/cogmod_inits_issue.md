# `cogmod_inits()` produces initial values Stan rejects in ~25% of chains

> **RESOLVED in cogmod 0.3.3** (2026-09-18, `dev` branch, commit `e920d44`).
> Both fixes recommended in section 7 were taken, and then some:
>
> - **The `A == 0` branch** (fix 1, the real bug). `cogmod_lognormal_acc_ltails()`
>   and `cogmod_lognormal_ldiff_Phi()` now go through one new function,
>   `cogmod_log_Phi()`: `erfc` in the body of the distribution and, below
>   `x = -25`, the asymptotic tail expansion, finite and differentiable as far
>   as `x = -1e150`. The difference of log-CDFs is taken with `log1m_exp()`
>   rather than as a quotient of two minute numbers. `_logcdf()` / `_logsurv()`
>   also compute only the tail they were asked for when `sigmabias = 0`, so a
>   discarded `-inf` node can no longer poison the tape. Measured: 6 of 72
>   gradients over a grid of decision times and sigmas were non-finite before
>   and none are now. The same route was applied to `cogmod_rdm()`, whose
>   gradient was finite but inexact (`std_normal_lcdf()`'s analytic partials).
> - **Tiered jitter** (fix 2). `z_*`, `zs_*`, `sd_*` and `sds_*` now get a
>   fifth of the population-level jitter (0.05 against 0.25), and `sds_*`
>   starts at 0.05 rather than the generic 0.25 — a smooth starts flat and
>   lets the data buy curvature.
> - **Fix 3** (`ndt` below every grouping level's own floor) was considered and
>   deliberately not taken; cogmod's `.ndt_start()` carries the reasoning, which
>   is that with the tiering in place the remaining exposure is under one trial
>   per data set.
>
> `fit_model.R` and `install_pkgs.R` both enforce cogmod >= 0.3.3. The fix is
> **not yet verified at production scale** — the failure rate grew with N, so
> the full-data run is the real test (`AGENT.md` 6.1, open question 7.2). The
> rest of this document is the original report, kept because it is the record
> of how the culprit was found and what the failure looked like.

Report for whoever works on `cogmod` (DominiqueMakowski/cogmod, installed
version 0.3.2, commit `e8eb18f6954e24d7edd8d42ca04414485181f8be`). Everything
below was measured on the Sussex Artemis cluster on 2026-09-17/18 with brms
2.21.0, cmdstanr 0.9.0, CmdStan 2.39.0, R 4.3.2 (GCC 12.3).

## 1. Symptom

A chain started from `init = cogmod_inits(f, data)` dies before its first
iteration with:

```
Chain 2 Rejecting initial value:
Chain 2   Gradient evaluated at the initial value is not finite.
Chain 2   Stan can't start sampling from this initial value.
Warning: Chain 2 finished unexpectedly!
```

The message is Stan's *gradient* check, not the log-density check, so the
log density at the init is finite and one or more partial derivatives are
`NaN`/`Inf`. Because brms passes explicit inits, CmdStan does **not** retry
with a fresh draw as it would for random inits: the chain is simply lost. When
every chain of a task is rejected, brms errors with
`Fitting failed. Unable to retrieve the metadata.` and the task produces
nothing.

## 2. How often, and where

Cold-start chains from the 2026-09-17/18 test matrix (`gam_lnr`, Muller-Lyer
data, the formula in section 5):

| participants | rows | chains started | rejected |
| --- | --- | --- | --- |
| 120 | 44,920 | 14 | 1 |
| 480 | 180,988 | 13 | 5 |
| 2,215 (all) | 323,981 | 2 | 1 |
| **total** | | **29** | **7 (24%)** |

Rejected chains: jobs 11368818 (480, 1 chain, whole task lost), 11371706 ch.2,
11371972 ch.2, 11371974 ch.1, 11371975 ch.3, 11372417 ch.2, 11373507 ch.2, plus
11368070_4 ch.1 in an earlier run. Logs are in
`/mnt/lustre/scratch/psych/dmm56/IGComputational/*.err`.

Two controls point at the jitter rather than at the target values:

- 12 chains started from `cogmod_inits(f, data, warmstart = <pilot table>)`,
  which uses `jitter = 0.05` around pilot posterior means: **0 rejected**.
- The rate does not fall with data size; if anything it rises (the number of
  participant-level random effects that get jittered grows with N).

`cogmod_inits()` defaults to `jitter = 0.25` (`.init_fun`, `.jitter_bounded`),
applied on the log scale for lower-bounded parameters
(`lower + (v - lower) * exp(rnorm(1, 0, 0.25))`) and additively for
unbounded ones.

## 3. What the init looks like for `cogmod_lnr()`

`.init_targets(cogmod_lnr())` returns
`mu = 0.7, nuone = 0.7, sigmazero = 0.5, sigmaone = 0.5, sigmabias = 0.5, ndt = 0.1, poutlier = 0.02`;
`ndt` is then replaced by `.ndt_start(Y)` = `0.5 * quantile(Y, 0.01)`.
Links (from the family object): `mu`, `nuone` identity; `sigmazero`,
`sigmaone`, `sigmabias` softplus; `ndt` log; `poutlier` logit.
`.init_value()` puts the linked target on `Intercept_<dpar>`, zeros on `b_*`,
and `.default_value()` on everything else: `lower + 0.25` for lower-bounded
scalars (so every `sd_*` and every smooth `sds_*` starts at 0.25) and 0 for
unbounded ones (all `z_*`, `zs_*`, `bs_*`).

Every one of those values is then jittered with `jitter = 0.25`. For the
`z_*` standardized random effects that means N(0, 0.25) starts for **every
participant on every distributional parameter**; for `Intercept_ndt` on the
log scale a factor `exp(N(0, 0.25))`, i.e. up to ~1.6x in either direction
at 2 SD.

## 4. Where a non-finite gradient can come from

The generated model block is (brms 2.21, cogmod stanvars):

```stan
sigmazero = log1p_exp(sigmazero);  sigmaone = log1p_exp(sigmaone);
ndt = exp(ndt);  poutlier = inv_logit(poutlier);
for (n in 1:N)
  target += cogmod_lnr_lpdf(Y[n] | mu[n], nuone[n], sigmazero[n], sigmaone[n], sigmabias, ndt[n], poutlier[n], dec[n]);
```

and `cogmod_lnr_lpdf` (cogmod stanvars) guards the obvious cases:

```stan
if (sigmazero <= 0 || sigmaone <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) return negative_infinity();
if (Y <= 0) return negative_infinity();
real t_adj = Y - ndt;
if (t_adj <= 0) return log(poutlier) + lp_out;      // RT below ndt is absorbed by the outlier component
...
return log_mix(poutlier, lp_out, lp_dec);
```

So `RT <= ndt` is **not** a hard boundary in the value: the density stays
finite through the outlier mixture. The gradient is another matter, and the
candidates are the branches that compute a quantity and then discard it:

- `cogmod_lognormal_ldiff_Phi` / the race integral: `u2 < u1 ? log(u1) + log1m(u2 / u1) : negative_infinity()`
  and `br <= 0 ? ... : negative_infinity()`. When `u1` or `u2` underflow to 0
  for an extreme `t_adj / sigma` combination, `log(u1)` and `u2/u1` evaluate
  `0/0` or `1/0` on the autodiff tape even though the ternary picks the
  constant branch; the value is fine, the adjoint is `NaN`.
- `log(log1p(A) / A)` and `log(series)` with `series <= 0` guarded only in
  the value.
- `log_mix(poutlier, lp_out, lp_dec)` with `lp_dec = -Inf`: value finite,
  but the derivative with respect to the inputs of `lp_dec` can be
  `0 * Inf`.

Section 6 shows which one fires: none of the guarded branches above, but
the *unguarded* `A == 0` shortcut in `cogmod_lognormal_acc_ltails`, which
calls Stan's `lognormal_lcdf` / `lognormal_lccdf` directly. The point for the
fix is that a jittered init only has to land **one** of ~300k trials in the
underflow region of those two functions.

## 5. Reproduction

Data: `data/illusion_part{1,2,3}.csv` from the repo, Muller-Lyer rows,
`Illusion_Difference <- abs(...)`, `Illusion_DifferenceZ` /
`Illusion_StrengthZ` as in `analysis/server/fit_model.R`. Formula:

```r
t2f <- function(lhs) as.formula(paste(lhs, "~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c('cr', 'cr')) + (1 | Participant)"))
f <- bf(RT | dec(Error) ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, k = c(5, 5), bs = c("cr", "cr")) + (1 | Participant),
        t2f("nuone"), t2f("sigmazero"), t2f("sigmaone"), sigmabias = 0, t2f("ndt"),
        poutlier ~ 1 + (1 | Participant), family = cogmod_lnr())
```

Any 480-participant subset reproduces it in roughly one chain out of two
with `init = cogmod_inits(f, data)`; 30 participants in about one draw out
of fifteen. Without sampling: compile the brms Stan code once with cmdstanr
(`make_stancode(...)` + `cmdstan_model()`), then
`mod$diagnose(data = standata, init = list(cogmod_inits(f, data)(i)), error = 1e6)`
for successive `i`. A rejected init makes `diagnose` fail with the same
"Gradient evaluated at the initial value is not finite" text (it does **not**
print the per-parameter table for a rejected init), which is what the
bisection over rows in `diag_bisect.R` works around. (cmdstanr's
`compile_model_methods = TRUE` route, which would give `grad_log_prob()`
directly, segfaults on this cluster's R 4.3.2 / GCC 12.3 combination.)

A unit-level reproduction needs no data at all: call
`cogmod_lnr_lpdf(3.89 | 0.005, 2.781, 0.239, 0.096, 0, 0.053, 0.023, 0)`
(first culprit) in a two-parameter Stan program and take its gradient;
`lognormal_lccdf(3.89 | -2.781, 0.096)` is `-Inf` with a `NaN` derivative.

Files on Artemis, all under
`/mnt/lustre/users/psych/dmm56/IGComputational/tests/inits_report/`:
`diag_bisect.R` + `bisect.out` (bisection, the culprit table above),
`diag_inits2.R` + `init_diagnose.csv` + `diag2.out` (rejection rates),
`stanvars.stan` (the LNR functions as generated by `cogmod_stanvars()`),
`lnr_model_plain.stan` (the full brms program).

## 6. Root cause: `cogmod_lognormal_acc_ltails()` at `A == 0` evaluates Stan's `lognormal_lcdf` / `lognormal_lccdf` in their underflow region

### 6.1 How it was found

Stan's `diagnose` method refuses a rejected init before printing anything, so
per-parameter gradients are not available for the failing draws. Instead the
init was held fixed and the **data** were bisected: the Stan data list was
sliced to subsets of rows (same participants, same smooth bases, so the
parameter vector is unchanged), and each subset was tested with
`mod$diagnose()`. Since the log density is a sum over trials, the subset that
still rejects contains the responsible trial. Every rejected init came down
to **one trial**, and that trial alone reproduces the rejection
(`single-row check: TRUE`). Script and log:
`tests/inits_report/diag_bisect.R`, `bisect.out`.

### 6.2 Rejection rates measured off-line (gradient at init, no sampling)

| participants | jitter | draws | rejected |
| --- | --- | --- | --- |
| 30 | 0.25 | 120 (two RNG streams) | 8 (6.7%) |
| 120 | 0.25 | 80 (two RNG streams) | 6 (7.5%) |
| 120 | **0.05** | 15 | **0** |
| 480 | 0.25 | 8 | 4 (50%) |
| production chains, 120-2215 | 0.25 | 29 | 7 (24%) |

### 6.3 The eight culprit trials

`z_los` is the losing accumulator's standardised log-time at the trial,
`(log(t_adj) - meanlog_los) / s_los` with `meanlog = -nu`, computed from the
init-implied distributional parameters at that row (intercept + smooth terms +
participant offset, on the link scale).

| subset | RT | dec | ndt (init) | t_adj | winner nu, s | loser nu, s | z_los | tail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 30 | 3.943 | 0 | 0.053 | 3.890 | 0.005, 0.239 | 2.781, 0.096 | **+43** | upper |
| 30 | 3.141 | 0 | 0.092 | 3.049 | 2.269, 0.244 | 1.974, 0.074 | **+42** | upper |
| 30 | 3.943 | 0 | 0.172 | 3.771 | 1.362, 0.824 | 1.912, 0.071 | **+46** | upper |
| 30 | 0.430 | 1 | 0.429 | 0.00099 | 0.866, 0.331 | 0.448, 0.138 | **-47** | lower |
| 30 | 0.474 | 1 | 0.465 | 0.00905 | 0.770, 0.352 | 0.111, 0.110 | **-42** | lower |
| 120 | 0.382 | 1 | 0.234 | 0.148 | 2.198, 0.827 | -0.618, 0.065 | **-39** | lower |
| 120 | 0.895 | 0 | 0.895 | 0.00027 | -1.261, 0.351 | -0.505, 0.153 | **-57** | lower |
| 120 | 0.443 | 0 | 0.438 | 0.0049 | 2.729, 1.171 | -0.048, 0.118 | **-45** | lower |

In every case the losing accumulator is 39-57 SD away from the observed
time. Nothing else is extreme: `poutlier` is 0.013-0.028, `t_adj > 0`, the
winner is at most 18 SD out (its `lognormal_lpdf` is finite there).

### 6.4 The mechanism

With `sigmabias = 0` (this project's formula, and cogmod's "plain LNR"),
`cogmod_lnr_decision_lpdf` calls `cogmod_lognormal_acc_ldens` and
`cogmod_lognormal_acc_logsurv`, and the latter goes through

```stan
vector cogmod_lognormal_acc_ltails(real t, real meanlog, real sigma, real A) {
  if (A == 0) {
    return [lognormal_lcdf(t | meanlog, sigma), lognormal_lccdf(t | meanlog, sigma)]';
  }
  ...  // A > 0: erfc-based tails, written precisely to avoid what follows
```

Stan's `lognormal_lcdf` / `lognormal_lccdf` are `log(0.5 * erfc(±z/sqrt(2)))`.
`erfc(x)` underflows to exactly 0 near `x = 26.5`, i.e. `|z| > ~37.5`, and
the value becomes `-Inf` while the partial derivative
`exp(-x^2) / erfc(x)` is `0/0 = NaN`. (The file's own header comment notes
that Stan's upper-tail log-CDF "is -Inf beyond y = 8.25" and switches to
`erfc` of a positive argument for the `A > 0` path; the `A == 0` shortcut
does not get that treatment.)

Two ways this reaches the gradient:

- **Upper tail (long RT, `z_los` = +42..+46):** `lognormal_lccdf` is the
  survival term that *is* used. `lp_dec = -Inf`, and
  `log_mix(poutlier, lp_out, -Inf)` gives a finite value
  (`log(poutlier) + lp_out`), so the log density passes Stan's check. Its
  gradient is `0 * NaN = NaN`.
- **Lower tail (fast RT or `t_adj` of a few ms, `z_los` = -39..-57):**
  `lognormal_lccdf` is fine (`log(1) = 0`) but the **discarded**
  `lognormal_lcdf` in the same two-vector underflows. Reverse-mode autodiff
  still runs the chain rule through the discarded node with adjoint 0, and
  `0 * NaN` poisons the adjoints of `t`, `meanlog` and `sigma`.

That is why the value is finite and only the gradient is not, why one trial
suffices, and why the rate grows with the number of trials.

### 6.5 Why the inits land there

`z_los` of 40+ needs a small `s_los` (0.065-0.15 here, against an intended
start of 0.5) together with a large gap between `log(t_adj)` and `-nu_los`.
Both come from the `jitter = 0.25` applied to the smooth coefficients
(`zs_*`, N(0, 0.25) each, multiplied by `sds_* = 0.25` and by tensor-product
basis values that reach tens) and to the participant offsets: at the culprit
rows the linear predictors have moved 1-2.5 units on the link scale
(`softplus^-1`: 0.5 -> 0.07; `log`: `ndt` 0.15 -> 0.89 s, above 98 of that
participant's 128 trials). With `jitter = 0.05` no init failed. Fewer
participants means fewer participant-specific `ndt`/`sigma` excursions and
fewer trials to hit, hence the size dependence.

## 7. Fixes, in order of payoff

1. **Fix the `A == 0` branch of `cogmod_lognormal_acc_ltails()` (the real
   bug).** Compute only the tail that is asked for, on the side where it is
   small, with a stable log-tail: `std_normal_lcdf(z)` for the lower tail and
   `std_normal_lcdf(-z)` for the upper tail are finite to `|z| ~ 37` and, past
   that, return `-Inf` with a finite derivative in recent Stan Math; or reuse
   the file's own `erfc`-of-positive-argument construction, clamped
   (`fmax(0.5 * erfc(x), 1e-300)`) so `log` never sees 0. Concretely:

   ```stan
   if (A == 0) {
     real z = (log(t) - meanlog) / sigma;
     return [std_normal_lcdf(z), std_normal_lcdf(-z)]';   // log F, log S
   }
   ```

   and, so that no `-Inf` with a `NaN` partial can survive on the tape, have
   `cogmod_lognormal_acc_logsurv()` / `_logcdf()` compute just their own
   element rather than both. This also removes the same hazard from every
   leapfrog step during warmup, where the sampler visits these tails
   routinely (the `NaN` gradient is what a divergence looks like there).
   `cogmod_rdm`, `cogmod_lba*` and the other race families that share this
   helper get the fix for free.
2. **Scale the jitter to the term.** 0.25 on `zs_*` smooth coefficients
   and on `z_*` participant offsets moves the linear predictor by 1-2.5 link
   units; 0.05 produced no rejection in 15 + 12 chains. A per-block jitter
   (e.g. 0.25 on intercepts and `b_*`, 0.05 on `zs_*` / `z_*` / `sds_*`)
   keeps between-chain dispersion where Rhat needs it and takes `ndt` and
   `sigma*` out of the corner. Independently of 1, this alone would have
   avoided all eight failures here.
3. **Start `ndt` below every participant's own floor.** `.ndt_start()` is
   `0.5 * quantile(Y, 0.01)` over all trials (0.15 s here); one culprit init
   had a participant at `ndt = 0.89 s`. Using the minimum over participants
   of each one's 1% quantile as the ceiling for the intercept-plus-offset,
   and jittering `Intercept_ndt` downward only, is cheap insurance.
4. **Make the init function self-checking** (optional, outside cogmod).
   `cogmod_inits()` returns a `function(chain_id)`; a wrapper that evaluates
   the gradient once (a compiled `cmdstan_model()` plus `$diagnose()`) and
   redraws on rejection would turn a lost chain into a 2-second retry. It
   needs the compiled model, which brms only exposes after fitting, so it is
   a documented pattern rather than a default.

## 8. Why it matters operationally

At full data a chain costs 12-20 h. Losing one chain in four is a quarter of
the cluster budget, and it forbids the otherwise attractive layout of one
chain x 16 threads per task (half the wall time per chain), because a single
rejection then loses the whole task. See `analysis/server/AGENT.md` 4.5 and
6.7.
