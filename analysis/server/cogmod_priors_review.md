# cogmod: default priors to change, from priors vs. posteriors on 18 full-data fits

2026-10-12. These are changes to propose for cogmod's default priors. They come
from a check of every combined fit of this project against the priors it was
sampled with. The next step of the project is to refit every model with longer
warmup and more iterations. Any prior change has to land in cogmod before that
refit, or it costs a second refit.

**Status 2026-10-12: cogmod 0.3.5 (`dev`) implements R1, R5 and R6; R2-R4 were
not implemented.** Before implementing, each proposal was refitted on cogmod's
vignette data on Artemis, current priors against proposed, with the same seeds
(`benchmarks/driftone_priors/README.md` in cogmod).

- **R1, R5, R6:** moved no posterior beyond noise and cost nothing.
- **R1 extended:** `cogmod_bisa()` got the same `boundary` prior. Its
  `boundary` was flat too, and missing from the table below.
- **R5:** the `bias` row also writes out `logistic(0, 1)`. Without it,
  `bias ~ 0 + Intercept` inherits the slope prior.
- **R2-R4 (`driftone`), not implemented:**
  - They predicted the same error rates and RTs.
  - They cost 9-26% more leapfrog steps per draw on the RDM.
  - On the LBA2 they raised the divergent chains from 1 in 12 to 5 in 12
    and 10 in 11.
  - The tail position they target is structural. Below a drift of a few
    hundredths the likelihood is flat, so a fenced plateau's posterior sits in
    its prior's tail, and no prediction depends on where. Expect it in
    `./hpc priors`.

The floor is raised to 0.3.5 here (`AGENT.md` §2).

**For the agent making the changes.** The cogmod repo is
`C:/Users/domma/Dropbox/Software/cogmod`, branch `dev`, at version 0.3.4
(`5382c79` at the time of writing). The cluster fits used `e8c0555` and
`732be30`. `R/core_choice.R` and `R/cogmod_priors.R` have not changed since
either. Items marked
*Recommended* should be done; the user has signed off on all six. The
*Decided against* and *Watch only* sections need no code. The *Project-side*
section is not for cogmod. Each change goes in `NEWS.md`, as for 0.3.4. The changes alter the
priors of fitted models, so they warrant a version bump (0.3.5); see
*After the changes*.

---

## How the check was done

`analysis/server/priors.R` runs on the cluster as `./hpc priors`, and the
results are shown in `analysis/1_models.qmd`. For every population-level
parameter of every combined fit, it records:

- the prior brms resolved for it, from `prior_summary(m, all = TRUE)` with
  brms's own inheritance (coefficient, then group, then class, within a dpar);
- the posterior median and 90% interval;
- the **prior quantile** of the posterior median: the prior's CDF at that
  median, so every parameter is on the same 0–1 scale;
- the **contraction**, `1 - var(posterior) / var(prior)`.

The flags are:

- **conflict**: the posterior median is beyond the prior's 1st or 99th
  percentile;
- **tail**: beyond the 5th or 95th percentile;
- **weak data**: contraction below 0.1;
- **flat prior**: brms's improper default.

**The resolved priors are checked against Stan's own `lprior`.** The summed
log prior densities reproduce the stored `lprior` at every draw, to 1e-5, in
all 18 fits. The 0.3.3 fits have one exception, a constant offset, which comes
from brms 2.21: it subtracts a truncated prior's normalising `lccdf` once per
*vector* rather than once per element. For example, a half-t on the three
`sds` of a `t2()` is off by `2 log 2`. The offset is constant, so it does not
affect sampling, and it is not a cogmod issue.

**The data.** About 324k trials per illusion: two-alternative choices with RTs
in seconds, error rate about 24% overall and 1–81% across conditions. There
are three illusions: Müller-Lyer, Vertical-Horizontal and Ebbinghaus. Every
free decision parameter gets
`t2(Illusion_Difference, Illusion_Strength, k = c(5, 5), bs = "cr")` plus
`(1 | Participant)`. `poutlier` gets `1 + (1 | Participant)`.

**The models.** The family is `cogmod_lnr()`, `cogmod_ddm()`, `cogmod_rdm()`
or `cogmod_lba2()`. LNR-6, DDM-5 and RDM-5 also free `sigmabias` or
`sigmadrift`. The Ebbinghaus fits, `gam_lnr_ebbinghaus` included, and both
`gam_lba` fits used cogmod 0.3.4. The others used 0.3.3 or older.

**Caveats.**
- **Reliable:** LNR, DDM-4, DDM-5 and RDM. Their R-hat is about 1.00–1.05.
- **Not converged:** LNR-6, RDM-5, LBA and the LNR Ebbinghaus refit. The first
  three have chains in different modes; LNR Ebbinghaus has chains that drift.
  Their pooled posteriors can be wider than the prior, so their flags are
  weak evidence and are kept apart below.
- **Large data:** with this much data, a conflict means the prior is pulling
  against the data, not that it decides the answer. The point of fixing these
  is a proper posterior, less pull, and defaults that do not fight a typical
  2AFC data set.

---

## Summary

| # | Change | Family | Status |
| --- | --- | --- | --- |
| R1 | Give `boundary` a prior (intercept and slopes); it is flat | DDM; same gap in `invgaussian`, `exwald` | **Recommended** |
| R2 | Re-centre and widen the `driftone` intercept prior | RDM | **Recommended** |
| R3 | Re-centre and widen the `driftone` intercept prior | LBA2 | **Recommended** (weaker evidence) |
| R4 | Widen the `driftone` slope prior | RDM, LBA2 | **Recommended** |
| R5 | Give the DDM `bias` slopes a prior; they are flat | DDM | **Recommended** |
| R6 | Move the `poutlier` intercept prior from `normal(-5, 1)` to `normal(-6, 1)` | all shifted families | **Recommended** (user's decision) |
| — | Scaled prior on the response's (`mu`) smooth slopes | all | **Decided against**: widened in the project instead |
| — | `sigmaone` / `sigmabias` intercepts | LBA2, LNR, RDM | Watch only, no change |
| — | `sds` scaling of 0.3.4 | all | **Confirmed**, keep |

---

## Recommended

### R1. `boundary` has no prior in `cogmod_ddm()`

**What happens.** `cogmod_priors()` leaves the DDM `boundary` intercept *and*
its blanket slope prior empty. brms then uses its flat default, an improper
prior. In all six DDM fits, 0.3.4 included, the Stan code has no `lprior` term
for `Intercept_boundary`. This project's `fit_model.R` fills the flat slopes
with a blanket normal (its patch, see *Project-side*), but not the flat
intercept.

**Where.** `R/core_choice.R`, the `cogmod_ddm` registry entry's `prior` list,
which names only `sigmadrift`, `sigmabias` and `sigmandt`. The link is
`softplus`; `init` has `boundary = 1`.

**The same gap exists elsewhere.** Checked with `cogmod_priors()` on the
`dev` source:

| family | `boundary` intercept | `boundary` slopes | `mu` slopes | other flat slopes |
| --- | --- | --- | --- | --- |
| `cogmod_ddm` | FLAT | FLAT | FLAT | `bias` |
| `cogmod_invgaussian` | FLAT | FLAT | FLAT | |
| `cogmod_exwald` | FLAT | FLAT | FLAT | |
| `cogmod_rdm`, `cogmod_lba2` | `normal(0, 1)` | `normal(0, 0.5)` | FLAT | |
| `cogmod_lba1` | `normal(0, 1)` | `normal(0, 0.5)` | FLAT | `sigma` |
| `cogmod_lognormal` | (none) | (none) | FLAT | `sigma` |

`mu` and `sigma` are the response's own parameters and are left flat on
purpose (see *Decided against*). `boundary` and `bias` are not response parameters, and
`?cogmod_priors` states the principle: "slope and group-level priors are
deliberately narrow" for the dpars a family names.

**Evidence.**

| fit | `boundary` intercept (softplus link) | natural scale |
| --- | --- | --- |
| DDM-4, three illusions | 1.37, 1.39, 1.41 | ≈ 1.6 |
| DDM-5, three illusions | 2.03, 2.11, 2.19 | ≈ 2.2–2.3 |

The posteriors are narrow (90% interval ± 0.03), so the likelihood identifies
the intercept and a proper prior will not move it. The boundary slopes are all
near 0 (|b| ≤ 0.24 under `normal(0, 1)`).

**Suggested.** Change the threshold-separation priors only where the boundary
means a separation; for the DDM, on the `softplus` link:

- intercept `normal(1, 1)`: median `softplus(1) = 1.31`, 95% between about
  0.32 and 3.0, which is the usual range for DDM separations at `s = 1`;
- slope `normal(0, 0.5)`, the slope prior the other threshold parameters
  already get;
- a natural-scale row for an omitted `boundary`, matched on the median as the
  file's convention requires (`?cogmod_priors`, *Parameters left out of the
  formula*), e.g. `lognormal(0.27, 0.5)`;
- a code comment, matching the existing ones, on why it is set: a flat prior
  on a `softplus` dpar is improper, and the separation trades off with drift
  and `ndt`.

`cogmod_invgaussian()` and `cogmod_exwald()` have the same flat `boundary`
(the threshold of a single-boundary Wald). It should get the same treatment,
but this project has no fit of either to calibrate against. Choose the centre
from those families' `init` and documentation, not from the numbers above.
Update the family list in `?cogmod_priors` ("A family may add rows of its
own... Seven do"), and `?cogmod_ddm` where it lists the priors.

### R5. DDM `bias` slopes have no prior

**What happens.** The `bias` intercept gets `logistic(0, 1)` (uniform on the
probability scale), but the blanket `b` for `bias` is empty, so its slopes are
flat. In this project they were filled with `normal(0, 1)` by `fit_model.R`.
The posteriors are within ±0.6 on the logit link, at prior quantiles of
0.36–0.72.

**Suggested.** Slope `normal(0, 0.5)` on the logit link, in the same registry
entry. The largest posterior here (0.58) sits at about the 88th percentile of
that, so it does not fight these data. Write it as R1's `boundary` row: the
registry entry needs a `bias` element with only a `slope` (and the
intercept's existing `logistic(0, 1)` kept), or the helper that turns the
registry into rows needs to allow that. Check how
`.priors_shifted()` builds rows before choosing.

### R2. RDM `driftone` intercept: `normal(3, 2)` sits far above every posterior

**What happens.** In all six RDM and RDM-5 fits the `driftone` intercept
posterior is at **-4.7 to -6.1** on its `softplus` link. That is at prior
quantiles of 3e-6 to 6e-5, the strongest conflict in the whole check. Five of
the six fits are clean (R-hat 1.00–1.04).

| fit | `driftone` intercept | natural (`softplus`) |
| --- | --- | --- |
| RDM Müller-Lyer | -6.07 [-7.71, -4.60] | 0.0023 |
| RDM Vertical-Horizontal | -6.07 [-7.17, -5.05] | 0.0023 |
| RDM Ebbinghaus (0.3.4) | -4.69 [-5.75, -3.85] | 0.009 |
| RDM-5 Müller-Lyer | -6.08 [-7.65, -4.59] | 0.0023 |
| RDM-5 Vertical-Horizontal | -5.72 [-6.89, -4.67] | 0.0033 |

The prior puts 99.9% of its mass above `softplus(-3.2) = 0.04`.

**Why the prior is off.** The registry comment says "a drift sits around 3
on this link". That is true of the *winning* accumulator: `mu`'s intercept is
1.6–2.1 here, a drift of about 1.8–2.2. But `driftone` is accumulator 1, and
with the usual `dec()` coding (1 = error) that is the error accumulator. At
the reference condition (centred predictors) its drift is near zero, and
errors there come from the diffusion noise of a nearly driftless accumulator.
The same comment notes that "a zero-drift accumulator still finishes", so this
is a legitimate region for an RDM, not a pathology.

**Note on contraction.** It is 0.78–0.92: the data inform the value, but the
prior still pulls it up by a noticeable amount.

**Suggested.** A prior that allows a near-zero drift without assuming it, for
example `normal(-2, 3)`:

- natural median `softplus(-2) = 0.13`;
- 95% from `softplus(-7.9) ≈ 4e-4` to `softplus(3.9) ≈ 3.9`.

That still covers the equal-drift case, and it is still proper and fences the
far tail. This data set's posteriors sit at z ≈ -1 to -1.4 under it. Also
revisit the natural-scale row `lognormal(1, 0.75)` (median 2.7) for the same
reason: it should allow drifts near zero, for example a lognormal with a
median of about 0.15 and an SD of about 1.5 on the log scale. Update the
registry comment, the paragraph on `cogmod_rdm()` in `?cogmod_priors`, and
the table in *Parameters left out of the formula*.

**Do not mirror this onto `mu`.** `mu` is the response's intercept and keeps
brms's default.

### R3. LBA2 `driftone` intercept: `normal(1, 2)` is the same mismatch

**What happens.** In both LBA fits (0.3.4) the `driftone` intercept is at
**-5.0** (Müller-Lyer) and **-5.9** (Vertical-Horizontal) on the identity
link. Those are prior quantiles of 0.001 and 0.0003. These fits are bimodal
(R-hat 1.65 and 2.20; `cogmod_lba_modes_issue.md`), but both modes are far
below the prior: the chains sit at about -4.5 and -6.8 in Müller-Lyer, and at
-5.5 and -9.8 in Vertical-Horizontal. So the conflict does not depend on which
mode is right.

**Why.** It is the same story as R2. The error accumulator's untruncated mean
drift is strongly negative, so that a positive drift is rare. In
`?cogmod_lba2` and in `cogmod_lba_priors_inits.md` (B2), `normal(1, 2)` was
described as "centred where LBA drifts sit on the `sigmazero = 1` scale".
That holds for the winning accumulator, not for the error one.

**Suggested.** `normal(-2, 3)`, as for the RDM. The identity link needs no
translation. It still fences the plateau the prior exists for: the docs say a
flat prior "lets the drift run off to -12 and beyond", and -12 is z = -3.3
under this prior. A side effect is relevant to the mode problem: the prior
penalty between the two Müller-Lyer modes drops from about 3.8 to about 0.9
log units. That is negligible against their 1,330 log-likelihood gap, so do
not expect it to fix the split. The natural-scale row is the same
`normal(1, 2)`, so change both.

*Weaker evidence than R2, because the fits are not converged. The conclusion
holds in every mode, so it is still recommended.*

### R4. `driftone` slope priors are too narrow for an error accumulator

**What happens.** The first `t2()` column of `driftone` conflicts, or nearly
does, in every RDM and LBA fit:

| fit | slope prior (after 0.3.4's rms scaling) | posterior | prior quantile |
| --- | --- | --- | --- |
| RDM Müller-Lyer (0.3.3) | `normal(0, 2)` | 4.65 | 0.990 |
| RDM Vertical-Horizontal (0.3.3) | `normal(0, 2)` | 5.88 | 0.998 |
| RDM Ebbinghaus (0.3.4) | `normal(0, 2.11)` | -4.14 | 0.025 |
| RDM-5 Müller-Lyer / Vertical-Horizontal (0.3.3) | `normal(0, 2)` | 4.74 / 5.50 | 0.991 / 0.997 |
| LBA Müller-Lyer (0.3.4) | `normal(0, 1.87)` | 4.39 | 0.991 |
| LBA Vertical-Horizontal (0.3.4) | `normal(0, 1.81)` | 5.44 | 0.999 |

The sign flips for Ebbinghaus because that illusion's strength runs the other
way. The second and third columns are fine.

**Why.** It is the same mechanism as R2 and R3. The error accumulator's drift
goes from about zero in the easy conditions to a real rate where the illusion
misleads. Error rates here run from 1% to 81% across cells. Near the
`softplus` floor that is a 5–10 unit swing in link space; on the LBA's
identity link it is a swing in the mean drift of a similar size. A slope prior
sized for "a condition effect worth a unit or two" (the registry comment)
fights it.

**Suggested.** `slope = "normal(0, 4)"` for `driftone` in both the
`cogmod_rdm` and `cogmod_lba2` entries. 0.3.4 then divides it by each
unpenalised column's rms. The largest posterior here would sit at about
z = 1.4. Update the registry comment ("a drift sits around 3...") in the same
edit as R2.

### R6. `poutlier` intercept: move the centre from -5 to -6

**Decision (user, 2026-10-12): change the default to `normal(-6, 1)`.**

**Evidence.** In all 18 fits the posterior is below the prior's centre, at
-6.0 to -8.2 on the logit link, i.e. 0.25% to 0.03% outliers:

- **Conflict in 9 fits:** all five RDM and RDM-5 fits, LNR Vertical-Horizontal
  and Ebbinghaus, LNR-6 Vertical-Horizontal and LBA Vertical-Horizontal.
  Prior quantiles are 0.0006–0.009.
- **Tail in 7 more.**
- **Unflagged:** only DDM-4 Müller-Lyer and Vertical-Horizontal (-6.0 and
  -6.2, quantiles 0.16 and 0.11).

The contraction is 0.95–0.99, so the data dominate, but the pull is
systematic and in the same direction everywhere.

**Why these data count as evidence for the default.** No early responses
were trimmed. `0_preprocessing.qmd` drops only the following:

- RT ≥ 6 s and RT ≤ 0. Its anticipation cut at 1/6 s is commented out.
- Whole blocks with a mean RT under 1/3 s or over 3 s, or with 45% errors or
  more.

So fast guesses are kept at the trial level, and the outlier component sees
them as they come. The estimates also agree with the reference the docs
already cite: Miller (2024) puts the outlier proportion *below 0.5%*, while
the current `-5` centre is 0.67%.

**What to change.**
- **Link scale:** `normal(-6, 1)` in `.SHIFTED_BASE_PRIORS`
  (`R/cogmod_priors.R`, line ~921). The centre is 0.25%, with about 95% of the
  mass between 0.035% and 1.7%. The most extreme posterior here (-8.2) would
  sit at z = -2.2.
- **Natural scale:** keep the `exponential` matched on the median, as the docs
  do now. That gives `exponential(280)`: median 0.0025, against
  `plogis(-6) = 0.0025`.
- **Text to update:**
  - the table and the paragraph in `?cogmod_priors` ("`normal(-5, 1)` is
    centred at about 0.7%...", and the *Parameters left out of the formula*
    table);
  - the comment in `R/cogmod_checkdata.R` (line ~343) that mentions the
    default;
  - any other place that quotes `normal(-5, 1)` or `exponential(100)`
    (grep for both).
- **Keep the *Trimmed data* advice** (`normal(-7, 0.5)`). It is still right
  for data that were trimmed.

---

## Decided against

### `mu`'s smooth slopes stay brms's (no change in cogmod)

cogmod leaves the response's slopes to brms, on purpose (`?cogmod_priors`:
"the response's slopes are left to `brms` because they are the effects being
estimated"). That includes the unpenalised `t2()` columns. This project
fills the flat `b` with a blanket prior in `fit_model.R`. That blanket,
`normal(0, 1)`, was too narrow for DDM-5: the first `mu` column is at
±2.2–2.5, prior quantiles 0.987, 0.992 and 0.005.

**Decision (user, 2026-10-12):** widen the project's blanket rather than add
smooth-specific rules to cogmod. Done in `fit_model.R` (see *Project-side*).
**Do not change cogmod for this.**

---

## Watch only: `sigmaone` and `sigmabias` intercepts

- LBA `sigmaone` intercept is 2.0–2.4 against `normal(0, 1)` (tail or
  conflict).
- `sigmabias` intercepts in LNR-6, RDM-5 and LBA are -2.2 to -3.3 against
  `normal(0, 1)`.

All of these come from non-converged, multimodal fits (R-hat 1.5–2.2), where
the pooled posterior mixes modes and is sometimes wider than the prior.
`sigmaone` also lives on the LBA's scale ray. **No change on this evidence.**
Revisit after the refit, when `./hpc priors` is re-run.

---

## Confirmed: keep as is

- **0.3.4's `sds` priors** (`exponential(rate)` per smooth, so `exponential(1)`
  in link units):
  - Under 0.3.3, the second smooth SD was in conflict across the LNR, DDM-5,
    RDM and LBA fits. Posteriors were 3–30 against `exponential(1)`, and 10–30
    against `mu`'s half-`student_t(3, 0, 2.5)`.
  - Under 0.3.4 (rate about 0.084–0.086 for these bases, mean about 12) there
    is no conflict left. Some smooth SDs, mostly `ndt`'s, sit in the low tail
    (0.2–0.5, prior quantiles 0.015–0.04). Their contraction is above 0.999,
    so those smooths are simply flat and the data say so.
- **0.3.4's rms-scaled slope priors** on the family-named dpars: fine apart
  from `driftone` (R4).
- **Participant SDs** (`exponential(1)`): all at prior quantiles 0.03–0.95.
  The only tail is the DDM `bias` SD of 0.09–0.16, under brms's
  half-`student_t(3, 0, 2.5)` in the 0.3.3 fits; it is harmless.
- **`ndt`** `normal(-1.2, 0.5)`, the LNR `nuone`, and the `sigmazero` /
  `sigmaone` / `sigmabias` rows of the converged models: no flags.

---

## Project-side (not cogmod)

This is listed so the cogmod agent does not try to do it.

1. **`mu` slopes: done 2026-10-12.** `analysis/server/fit_model.R` now fills
   flat blanket `b` priors with `normal(0, 2)` instead of `normal(0, 1)`.
   Once R1 and R5 land, `boundary` and `bias` get their slope priors from
   cogmod, and the loop skips them by itself, so the blanket reaches only
   `mu`.

---

## After the changes

1. Add a `NEWS.md` entry for each change, saying which models' priors change.
   Bump to 0.3.5. The project's `fit_model.R` and `install_pkgs.R` enforce a
   floor (`IGC_COGMOD_MIN`, now 0.3.4). Raise it so the refit cannot run on
   the old priors (`AGENT.md` §2 describes the mechanism).
2. Check with `cogmod_priors()` on a formula that gives every dpar a smooth:
   no `Intercept` or blanket `b` row of a family-named dpar should be empty,
   for the DDM, invgaussian, exwald, RDM and LBA2. The script used for the
   table in R1 is a loop of `cogmod_priors(bf(..., family = f()), d)` over the
   families, listing empty `prior` rows by class.
3. Run cogmod's tests.
4. Project side: `./hpc install cogmod`, but never while an array has pending
   tasks (`AGENT.md` §2). Then the refit. After it, `./hpc priors` and
   `./hpc pull priors`, and re-render `1_models.qmd` to confirm the flags are
   gone.
