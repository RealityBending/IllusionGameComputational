# IllusionGameComputational

Bayesian cognitive models (LNR, DDM) of reaction times and errors in a visual
illusion task, fitted with brms + cmdstanr and the
[cogmod](https://github.com/DominiqueMakowski/cogmod) package.

## Model fitting runs on a cluster, not here

Fits take 12-20 hours per chain on 324k rows and are run on **Artemis**, the
University of Sussex HPC. Do not try to fit them locally.

Cluster-level instructions (access, storage, partitions and quotas, the R/Stan
toolchain, job conventions, troubleshooting, housekeeping) live in the **lab
HPC hub**: <https://github.com/RealityBending/Lab/tree/main/hpc> — start at
`hpc/README.md`, and follow its rules for agents. Everything specific to this
project lives in `analysis/server/`:

- **`analysis/server/README.md`** — the command reference. Read it before
  running anything against the cluster.
- **`analysis/server/AGENT.md`** — the measurements, the traps, and the
  reasoning behind the current settings. **Read this before changing a
  partition, a resource request, warmup, or the number of chains.** Several of
  those settings look arbitrary and are not; each one that was changed
  carelessly cost a failed run, and they are written up with the evidence.
- cogmod's [`NEWS.md`](https://github.com/DominiqueMakowski/cogmod/blob/main/NEWS.md)
  (0.3.3 entries) records the cold-start initialisation failures and their fix.
  It also explains why freeing `sigmabias` / `sigmandt` makes the DDM so much
  dearer per gradient: a branch into numerical quadrature, not model geometry.
  Read that, and `analysis/server/AGENT.md` §4.7.1, before proposing anything
  about the seven-parameter DDM.
- `analysis/server/cogmod_lba_modes_issue.md` — why the full-data `gam_lba`
  fit split into two posterior modes: not a code bug, but the error
  accumulator's identifiability where it rarely wins. Also what cogmod could
  change, and the reparametrised variants. Read it before refitting any
  `gam_lba*`.
- `analysis/server/cogmod_lba_priors_inits.md` — the cogmod defaults that were
  to be fixed before re-running `gam_lba`, most of them now changed in cogmod
  0.3.4 on `dev` (its status table), and the full `gam_lba` specification.

The workflow, in one line:

```bash
cd analysis/server && ./hpc push && ./hpc fit <model>
```

`./hpc` with no arguments prints its commands; `./hpc models` lists the models.

## Ground rules for this repo

- **One job fits one model.** `analysis/server/models.R` is the single
  definition of what a model is — adding one is a single entry there and
  nothing else. Do not hard-code a formula into a fitting script.
- **Not every model in the registry should be submitted.** `gam_ddm7` is known
  not to be viable **at full data** (README → "`gam_ddm7`: do not submit it at full data",
  and `AGENT.md` §4.7). It is viable on a subsample — roughly 200 participants
  — and that run is worth doing; what is ruled out is the full 2,215. Check the
  README's "Who runs what" table before launching anything.
- **Never hard-code an account or a path.** Everything is an `IGC_*` variable
  with a default (see README → Paths). A second cluster account sets its own in
  `analysis/server/hpc.local`, which is gitignored.
- **Test runs get their own `IGC_MODELS_DIR`.** Fits use
  `file_refit = "never"`, so a leftover 30-participant test shard in the
  production directory is silently adopted by a production run.
- **The GlobalProtect VPN must be connected** for anything touching the
  cluster; the other cluster-wide rules (few SSH connections, nothing deleted
  without approval, ...) are in the hub's `README.md`.
- `analysis/models/` is gitignored; fitted `.rds` files are not committed.
- **`1_modelcomparison.qmd` does not load fits.** It plots the prediction files
  that `./hpc predict <model>` writes (README → "Predictions for the model
  comparison"). Anything it computes from a fit belongs in
  `analysis/server/predictions.R`, which the qmd and the cluster job share;
  the qmd itself only relabels and plots. Bump `igc_predictions_version` when
  the structure of what that file returns changes.
- **Never run brms post-processing on a cluster fit with the laptop's brms.**
  The fits are brms 2.21.0, and the laptop's 2.23.1 rebuilds their smooth
  bases with the wrong signs. `log_lik()`, `loo()`, `posterior_*()`,
  `fitted()`/`predict()` and `modelbased::estimate_*()` then return wrong
  numbers without any warning. Predict on the cluster, use brms 2.21.0, or
  call `analysis/server/keep_stored_basis.R` on the fit straight
  after `readRDS()`. `standata(m)` is unaffected, so it is not a check.
  `analysis/server/AGENT.md` §3.9 has the details.
- **New fits use cogmod 0.3.4** (decision 2026-10-02), which changes the
  priors of every model with a smooth. Whether to refit the older fits is
  decided once every wanted model is fitted, so do not propose it before then.
  The guards are a 0.3.4 floor, a per-shard `m$cogmod` stamp, and a combine
  step that refuses mixed versions. Read `analysis/server/AGENT.md` §2 before
  `./hpc install cogmod`.

## Two ways to fit more than one model at a time

**`sussexneuro`.** `dmm56` belongs to the `artemis_sussexneuro` group
(verified 2026-09-22), and that departmental partition carries its own
256-CPU group quota, separate from `long`'s per-user 140. A job there runs
*in parallel* with the production arrays at no cost to them —
`./hpc fit <model> --partition=sussexneuro`. Read the hub's
`artemis.md#sussexneuro` first (a shared group pool, over-quota jobs rejected
rather than queued, 60 days is not 60 usable days), and `AGENT.md` §4.2.1 for
how this project uses it.

**A second account.** The per-user CPU quota is the binding constraint, and it
is per account, so two people can fit different models concurrently. Setup is
one file and one line — see README → "Running from a second cluster account".
The model registry is shared through git; only the account is local. A model
one person fits must still be defined and committed in `models.R`.
