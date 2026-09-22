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
- `analysis/server/cogmod_inits_issue.md` — the cold-start initialisation
  failures and their root cause, fixed in cogmod 0.3.3.
- `analysis/server/cogmod_ddm_cost_issue.md` — why `gam_ddm7` costs 55x per
  gradient (a branch into numerical quadrature, not model geometry), what does
  and does not fix it, and the change cogmod would need. Read it before
  proposing anything about the seven-parameter DDM.

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
