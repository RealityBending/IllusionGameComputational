# `cogmod_ddm()` costs 18-55x per gradient when a between-trial variability is freed

> **Status 2026-09-22.** The cost is confirmed and its cause is understood and
> verified against source; **no fix has been made**. What follows is a report
> for whoever works on `cogmod` (DominiqueMakowski/cogmod), plus the record of
> what the cost rules out for this project.
>
> Two provenance notes, so nothing here is taken for more than it is:
>
> - **The benchmark table in section 2 was measured on 2026-09-18** and is
>   reproduced from `AGENT.md` §4.7.1, which is where it was first written
>   down. **The harness that produced it was never committed and is not
>   recoverable** — it is in no commit, no stash and neither dangling commit in
>   this repository. Section 5 specifies a harness that would re-derive the
>   same quantity; it is new, written 2026-09-22, and **has not been run**, so
>   it is not the provenance of the numbers above it and must not be presented
>   as such.
> - **Everything in sections 3 and 4 was read off cogmod source on 2026-09-22**
>   (version 0.3.3, branch `dev`, commit `9f9f809`) rather than inferred from
>   behaviour, and the line references are to that commit.

Measured against brms 2.21.0, cmdstanr 0.9.0, CmdStan 2.39.0, R 4.3.2
(GCC 12.3) on the Sussex Artemis cluster.

## 1. Symptom

`gam_ddm7` — the seven-parameter DDM, all three between-trial variabilities
free and each carrying a 2-D smooth — did not reach iteration 100 in 50 minutes
at **30 participants / 3,806 rows**, where every other model in the registry
finished all 400 iterations in 6-11 minutes (`AGENT.md` §4.7). That is upwards
of 25x the per-iteration cost on 1/74th of the production data.

The first guess was geometry: three weakly-identified variance parameters, each
with 25 tensor coefficients and a participant intercept, would plausibly make
the sampler take maximum-length trajectories and still not move. **That guess
was wrong**, and the distinction matters operationally: a geometry problem is
fixed by reparameterising and might be bought off with more hours, whereas this
one cannot be bought off at all.

## 2. The measurement

Measured 2026-09-18 with a standalone Stan program calling `wiener_lpdf` in
each of its forms, cost normalised by `n_leapfrog__` so the variants are
comparable.

### 2.1 The benchmark

| density variant | µs/obs/gradient | x the DDM-5 path |
| --- | --- | --- |
| classic 4-parameter (`gam_ddm4`) | 1.97 | 0.36 |
| 5-parameter, `sigmadrift` free (`gam_ddm5`) | 5.56 | 1.0 |
| one of `sigmabias`/`sigmandt` free (1-D quadrature) | 100.0 | **18.0** |
| both free (2-D quadrature, `gam_ddm7`) | 307.1 | **55.3** |
| both free but driven to `1e-5` | 144.7 | 26.0 |

The three tiers are not an artifact of the benchmark: they correspond exactly
to the three branches the source actually has (3.1).

### 2.2 It reproduces the fits it is explaining

The benchmark predicts `gam_ddm4` at 5.2 min against the 7.1 measured, and
`gam_ddm7` reaching ~9 iterations in 46 min against the observed "did not reach
100". At full data with 8 threads the classic path predicts 80 ms per gradient
against the 78 ms §4.4.1 measured for `gam_lnr`.

## 3. Root cause: an exact-zero branch into Stan's adaptive quadrature

### 3.1 The branch, as cogmod 0.3.3 actually writes it

`cogmod_ddm_decision_lpdf()` (`R/model_ddm.R:1437`). Abridged, but the control
flow and the tests are verbatim:

```stan
real sw = sigmabias * fmin(2 * w, 2 * (1 - w));

if (sw == 0 && sigmandt == 0) {
  if (sigmadrift == 0 && t > 0.01 * square(boundary)
      && cogmod_ddm_log_density_scale(t, v, boundary, w, 0, 0, 0) > -600) {
    return wiener_lpdf(y | boundary, tau0, w, v);              // wiener4
  }
  return wiener_lpdf(y | boundary, tau0, w, v, sigmadrift);    // wiener5
}
if (cogmod_ddm_log_density_scale(t, v, boundary, w, sigmadrift, sw,
                                 sigmandt) < -600) {
  return negative_infinity();
}
return wiener_lpdf(y | boundary, tau0, w, v, sigmadrift, sw, sigmandt,
                   1e-3);                                      // wiener7
```

Three paths, two **nested exact-zero tests**, and the benchmark's three tiers
line up with them one for one:

| branch | entered when | benchmark row |
| --- | --- | --- |
| `wiener4` | `sigmabias`, `sigmandt`, `sigmadrift` all exactly 0 | 1.97 |
| `wiener5` | `sigmabias`, `sigmandt` exactly 0; `sigmadrift` free | 5.56 |
| `wiener7` | either of `sigmabias`, `sigmandt` non-zero | 100.0 / 307.1 |

`sigmadrift` is marginalised analytically inside `wiener5` — a Gaussian
integral — which is why `gam_ddm5` is cheap. The other two have no closed form
and `wiener7` delegates them to **Stan's adaptive numerical quadrature**;
cogmod's own comment above that line says so. The `1e-3` is
`.DDM_WIENER_PRECISION` (`R/model_ddm.R:930`), the tolerance handed to the
adaptive rule.

Note that `sw == 0` is equivalent to `sigmabias == 0`, since
`fmin(2w, 2(1-w)) > 0` for any `w` strictly inside `(0, 1)`.

### 3.2 Why tight priors cannot reach the fast path

They cannot, and this is the part most likely to be re-discovered the expensive
way. The test is for **exact** zero; an estimated parameter never satisfies it,
whatever its prior. Driving both variabilities to `1e-5` still costs 26x — and
`1e-3` and `1e-5` measure the same, because both take the identical branch and
merely make the integrand narrow.

There is a second trap layered on the first: **"zero" is on the link scale.**
`sigmabias` is logit-linked, so η = 0 is *half* the maximum start-point range;
`sigmandt` is log-linked, so η = 0 is **one second** of non-decision-time
range. A prior centred on zero there is pathological rather than conservative.
cogmod's own `normal(-2, 1)` and `normal(-3, 1)` are already the sensible
locations.

### 3.3 Freeing one costs 18x, either one

Freeing exactly one variability costs 18.0x whichever one it is — the two
measure the same to three figures, because both enter the same 1-D quadrature.
So **which one to free is a modelling choice, not a performance one.**

If one is ever freed it should be **`sigmabias`**. Freeing `sigmandt` silently
redefines `ndt` as the *lower bound* of a between-trial distribution rather
than its midpoint (cogmod documents this), and `ndt` carries its own 2-D smooth
in every model in the registry — so the `ndt` surface would quietly stop being
comparable across the model set.

## 4. Fixes

### 4.1 The real fix: a fixed-node Gauss-Legendre rule, and its design constraints

Replace the delegation to Stan's adaptive quadrature with cogmod's own
fixed-node Gauss-Legendre rule. **The pattern already exists in the package**
and is the template to copy, not a new design:

- `.GAUSS_LEGENDRE` (`R/core_shifted.R:1426`) builds 64 nodes and weights once,
  by Golub-Welsch on the Jacobi matrix.
- `.pwald_sv()` (`R/core_shifted.R:1383`) is the R side: it marginalises the
  fixed-drift Wald over a zero-truncated normal drift on
  `[max(mu - 10 sigma, 0), mu + 10 sigma]`.
- `cogmod_wald_sv_lquad()` (`R/core_shifted.R:1581`) is the Stan side, and its
  nodes and **log-weights are pasted in as literals generated from the same
  `.GAUSS_LEGENDRE`**, so R and Stan integrate over literally the same points.
  Terms are assembled with `log_sum_exp` so the survival keeps its digits where
  every term is tiny.

The design constraints such a rule has to meet for the DDM:

1. **Two dimensions, not one.** `sigmabias` marginalises over the start point
   and `sigmandt` over the non-decision time, so the seven-parameter form needs
   a product rule. The `wiener7` row (307.1) is roughly 3x the 1-D row (100.0),
   which is the signature of a nested adaptive rule and sets the bar to beat.
2. **Node count is the whole question, and the Wald precedent argues both
   ways.** §4.7.1 projects a 3x3 rule at ~6x faster. The Wald's own note is
   that 64 nodes reach ~5e-15 and that **48 already suffice everywhere but the
   widest** configurations — in *one* dimension, over a *truncated normal*.
   That looks discouraging for 3x3, but the DDM's two integrands are
   **uniform** (a uniform start-point range and a uniform non-decision-time
   range), not truncated normal, and low-order Gauss-Legendre on a smooth
   function against a uniform weight is exactly the friendly case. The honest
   statement is that 3x3 is plausible and unverified; the node count must be
   chosen by measuring against the adaptive rule across the parameter box the
   fits actually occupy, as `.pwald_sv()`'s note documents having done.
3. **Gradients, not just values.** The reason cogmod routes every normal tail
   through `cogmod_log_Phi()` rather than `std_normal_lcdf()` is that the
   built-in's *value* is right and its *partials* are not. A fixed-node rule is
   differentiated node by node, which is precisely what makes it well-behaved
   here — but the same scrutiny applies: check the gradient, not only the
   density.
4. **It must not move the existing branches.** `wiener4` and `wiener5` are
   correct and cheap; only the `sw != 0 || sigmandt != 0` fall-through should
   change. The `-600` log-density guard in front of it exists because the
   seven-parameter form has the same `-inf`-with-`NaN`-partials failure as the
   classic one, and must survive whatever replaces the call.
5. **The payoff is what makes a six-parameter DDM routine.** A ~6x reduction is
   what would bring a six-parameter DDM into range as an ordinary fit, and a
   seven-parameter one into range in combination with a warm start.

### 4.2 A separate matter: the `sds` prior on a smoothed dpar

Recorded here because `AGENT.md` §5.7 cites this section for it.

`sds` is the smooth-wiggliness SD. On a logit- or log-linked dpar, brms's
default `student_t(3, 0, 2.5)` — a half-t with median ~1.9 on the link scale —
lets the smooth alone walk a `sigmabias` across its whole range or move a
`sigmandt` by a factor of seven, undoing the tight intercept the family put
there on purpose. It is the loosest prior in these models.

**This is fixed in cogmod 0.3.3, and the fix is subtler than "set a prior".**
`R/cogmod_priors.R:920-933` documents the mechanism: a smooth's `sds` arrives
the opposite way round from a group-level `sd` — brms fills the *blanket* row
itself and leaves the per-term rows empty, so the empty rows were "covered" and
the filled one was never a candidate. The result was that **until 0.3.3 an
`sds` on a dpar silently kept brms's default while `?cogmod_priors` said
`exponential(1)`.** 0.3.3 takes the blanket row and replaces it, so `sds` now
genuinely gets `exponential(1)`. The response's own smooth (`dpar == ""`) is
deliberately left alone, for the same reason the response's slopes are.

Consequence for this project: **whether a given fit has the tight or the loose
`sds` depends on which cogmod built its Stan code**, so it is a property of the
cluster library at submission time, not of anything in this repo. Check before
comparing surfaces across models fitted at different times.

### 4.3 What does not work, so no one spends a run finding out

- **Tight priors near zero.** 55x → 26x and no further (3.2).
- **Making the variabilities intercept-only.** The cost is per-observation and
  independent of how many coefficients feed the linear predictor, so dropping
  the smooths on the three variabilities changes nothing.
- **A warm start.** It buys back warmup, not per-gradient cost; the 500
  retained draws alone are ~18 days at full data.
- **More hours.** See section 6.

## 5. Re-deriving the benchmark

**Not the provenance of section 2.1** — the original harness was not preserved.
This is a specification for reconstructing the measurement, written 2026-09-22
and **not yet run**. Anyone acting on it should treat its output as a fresh
measurement and say so.

What has to hold for the result to be comparable with 2.1:

1. **Normalise by `n_leapfrog__`, not by wall time or by iteration.** The
   variants differ in how far the sampler moves per iteration, so only per
   gradient is comparable. Sum `n_leapfrog__` from the sampler diagnostics and
   divide wall time by it and by N.
2. **Drive the branch from data, not from a Stan flag.** Each variant must be a
   separate compiled program in which the relevant argument is a literal `0` or
   a parameter, because the branch under test is an exact-zero test at run
   time — passing a parameter that happens to equal zero takes the *fast* path
   and measures nothing.
3. **Five programs**, matching 2.1's rows: `wiener4`; `wiener5`; one of
   `sigmabias`/`sigmandt` free; both free; both free with the estimated values
   pinned near `1e-5`.
4. **One observation set, reused.** Same simulated RTs, same `boundary`, `w`,
   `v`, same seed across all five, so the only difference is the branch.
5. **Match the tolerance.** Pass `.DDM_WIENER_PRECISION` (`1e-3`) to the
   seven-parameter call, or the adaptive rule is being measured at a tolerance
   the package does not use.
6. **Report the parameter box.** The adaptive rule's cost depends on the
   integrand's width, so a single number is only meaningful with the
   `sigmabias`/`sigmandt` values it was measured at.

The natural place to run it is `short` at 30 participants; it is minutes of
compute, not hours.

## 6. What it means operationally

**Subsampling is the only thing that works.** The cost is per row, so per chain
at warmup 1000 + 500 draws, 16 threads, 1 chain per task:

| | full data | fits ~5 days on `long` |
| --- | --- | --- |
| `gam_ddm7` as specified | ~54 days | ~240 participants |
| `gam_ddm7`, the two pinned tight | ~26 days | ~500 participants |
| a six-parameter DDM | ~18 days | ~700 participants |

At 200 participants `gam_ddm7` is 2.7 days per chain if leapfrog settles at
255, and 4.7 if at 511 — against 8.7, over `long`'s wall, if treedepth stays
pinned at 10. **That last case is itself the finding**, so the run is
informative whether or not it completes. Give it **1 chain x 16 threads**
(§4.4.1: threading scales, more chains per task does not) and **its own
`IGC_MODELS_DIR`**.

The `sussexneuro` partition's 60-day ceiling does **not** change this. It
removes the wall, but 54 days of sampling for one model is not a good use of a
quota shared with the rest of the group, and the 200-participant run answers
the same question in under a week. See `README.md` → "`sussexneuro` — a second,
separate allowance".

`gam_ddm4` and `gam_ddm5` are unaffected: both take branches that never reach
the quadrature.

---

**Related.** `AGENT.md` §4.7.1 (the same finding in its operational context and
where the benchmark was first recorded), §4.7 (the smoke tests this explains),
§5.7 (priors); `README.md` → "`gam_ddm7`: do not submit it at full data";
`cogmod_inits_issue.md` (the other cogmod report, resolved in 0.3.3).
