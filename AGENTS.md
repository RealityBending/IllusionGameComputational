# IllusionGameComputational

Bayesian cognitive models (LNR, DDM) of reaction times and errors in a visual
illusion task, fitted with brms + cmdstanr and the
[cogmod](https://github.com/DominiqueMakowski/cogmod) package.

## Model fitting runs on a cluster, not here

Fits take 12-20 hours per chain on 324k rows and are run on **Artemis**, the
University of Sussex HPC. Do not try to fit them locally.

Everything for that lives in `analysis/server/`:

- **`analysis/server/README.md`** — the command reference. Read it before
  running anything against the cluster.
- **`analysis/server/AGENT.md`** — the measurements, the traps, and the
  reasoning behind the current settings. **Read this before changing a
  partition, a resource request, warmup, or the number of chains.** Several of
  those settings look arbitrary and are not; each one that was changed
  carelessly cost a failed run, and they are written up with the evidence.
- `analysis/server/cogmod_inits_issue.md` — the cold-start initialisation
  failures and their root cause, fixed in cogmod 0.3.3.

The workflow, in one line:

```bash
cd analysis/server && ./hpc push && ./hpc fit <model>
```

`./hpc` with no arguments prints its commands; `./hpc models` lists the models.

## Ground rules for this repo

- **One job fits one model.** `analysis/server/models.R` is the single
  definition of what a model is — adding one is a single entry there and
  nothing else. Do not hard-code a formula into a fitting script.
- **Never hard-code an account or a path.** Everything is an `IGC_*` variable
  with a default (see README → Paths). A second cluster account sets its own in
  `analysis/server/hpc.local`, which is gitignored.
- **Test runs get their own `IGC_MODELS_DIR`.** Fits use
  `file_refit = "never"`, so a leftover 30-participant test shard in the
  production directory is silently adopted by a production run.
- **The GlobalProtect VPN must be connected** for anything touching the
  cluster. `kex_exchange_identification: Connection reset` means sshd is
  rate-limiting a burst of connections, not that the VPN dropped — wait a
  couple of minutes rather than retrying in a loop.
- `analysis/models/` is gitignored; fitted `.rds` files are not committed.

## Working from a second cluster account

The per-user CPU quota is the binding constraint, and it is per account, so two
people can fit different models concurrently. Setup is one file and one line —
see README → "Running from a second cluster account". The model registry is
shared through git; only the account is local. A model one person fits must
still be defined and committed in `models.R`.
