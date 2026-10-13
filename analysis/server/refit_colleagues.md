# Refitting your models on your own Artemis account (October 2026)

Every model is being refitted under a new setup, and you are taking some of the
slowest ones on your own account. These need little attention: about four
check-ins over roughly two weeks.

What changed since the last fits:

- **cogmod 0.3.5** (`dev`), with new priors (`poutlier`, the DDM thresholds),
  LBA density fixes and new LNR starting values.
- **New data**: blocks with 1/3 or more errors on facilitating trials are
  excluded. The fits read it straight from GitHub, so you have nothing to do.
- **Longer runs**: 1,500 warmup + 1,000 draws per chain, 8 chains per model.
  These are the defaults in `fit_model.R` now, so do not set `IGC_WARMUP` or
  `IGC_SAMPLES`.

All commands run from `analysis/server/` in your clone of the repo, with the
**GlobalProtect VPN connected**. `./hpc` with no arguments lists the commands,
and `README.md` in the same folder is the full reference.

## Who fits what

| | models, in submission order | per shard on a zen5 node |
| --- | --- | --- |
| A | `gam_lba`, `gam_lba_ebbinghaus`, `gam_lba_verticalhorizontal`, `gam_rdm5_ebbinghaus` | LBA ~4.3 days (up to ~5.1), RDM-5 ~4.3 (up to ~5.7) |
| B | `gam_lnr6`, `gam_lnr6_ebbinghaus`, `gam_rdm5_verticalhorizontal`, `gam_lnr6_verticalhorizontal` | LNR-6 3.5-5.5 days (up to ~6.6), RDM-5 ~3.9 (up to ~4.8) |

Dom's message tells you which row is yours. Each model is 4 shards (array
tasks) of 2 chains on 16 CPUs. Your account can run 8 tasks at once on
`long` (its 140-CPU cap at 16 CPUs a task), i.e. two models. So the first two
models run together, and the other two start as their shards finish. That
makes about **10-12 days** for all four, or longer if the zen5 nodes are busy.
Submit all four at once anyway: waiting tasks start by themselves and cost
nothing.

## 1. Update (5 minutes)

```bash
git pull
cd analysis/server
./hpc check            # VPN + SSH
./hpc queue            # must list no igc_ jobs before step 2
```

**First time on this project?** Follow README → "Running from a second cluster
account" first (one `hpc.local` line, then `./hpc setup`, `./hpc install`,
`./hpc precompile`). If you have fitted a model here before, skip it.

## 2. Move your old fits out of the way

```bash
./hpc archive
```

This moves everything in your `models/` (shards, `combined/`,
`predictions/`, `priors/`) into `models/old_<today's date>/`. **Nothing is
deleted.** It matters because a fit reuses a shard of the same name if one is
already on disk. An old shard would now stop the job ("is from another run")
rather than be reused, but that is still a failed job. `archive` refuses to
run while any `igc_` job of yours is queued or running.

## 3. Install cogmod 0.3.5

```bash
./hpc push             # copies the new scripts to the cluster
./hpc install cogmod   # ~5 min on a compute node
```

The output must end with `cogmod      OK 0.3.5` and `cogmod >= 0.3.5 OK`. If
you have not run a plain `./hpc install` since late September, run it once
too: the predictions in step 7 need `modelbased`.

Your R library is shared by every project on your account, so do this only
when **no job of yours, from any project, is pending**. Once your fits are
submitted, **do not reinstall cogmod** until all of them are done. A task
that starts after a reinstall would be fitted on different priors from its
siblings, and the combine step would refuse it.

## 4. Smoke test (10-15 minutes)

Use a 30-participant run on the 2-hour `short` partition, in its own
directory, so it can never be mixed into the real fits. Replace `<you>` with
your cluster user name (and `psych` with your group, if that differs).

```bash
IGC_MODELS_DIR=/mnt/lustre/users/psych/<you>/IGComputational/smoke_refit \
IGC_NPARTICIPANTS=30 IGC_WARMUP=300 IGC_SAMPLES=100 \
  ./hpc fit gam_lba,gam_rdm5 --array=1 --partition=short --cpus-per-task=8 --mem=16G
```

B: the same with `gam_lnr6,gam_rdm5`. Then, with each job number `sbatch`
printed:

```bash
./hpc log <job number>
```

You want `participants: 30  rows: 3734`, `SUCCESSFUL`, a `REPORT` line and
`status 0`. On Dom's account these took 4-17 minutes, LNR-6 the longest. A
compile error within the first minute is the brms problem in `AGENT.md` §3.8.
Tell Dom rather than working around it.

## 5. Submit the real fits

A:

```bash
./hpc fit gam_lba,gam_lba_ebbinghaus,gam_lba_verticalhorizontal,gam_rdm5_ebbinghaus --constraint=amd_zen5
```

B:

```bash
./hpc fit gam_lnr6,gam_lnr6_ebbinghaus,gam_rdm5_verticalhorizontal,gam_lnr6_verticalhorizontal --constraint=amd_zen5
```

The order is the order they start in: the slowest first.

Use exactly that: no `--time`, no other partition, and keep
`--constraint=amd_zen5`. The zen5 nodes are 2-4 times faster than the
others, and these models only fit inside `long`'s 8-day limit on them. A
shard that hits the limit is lost.

**Note the four job numbers** that `sbatch` prints.

`./hpc queue` should then show 16 tasks. Some will be `PENDING`, which is
normal:

- `QOSMaxCpuPerUserLimit`: waiting for your own earlier tasks to finish
  (the last two models).
- `Resources` or `Priority`: waiting for a zen5 node. There are only three,
  shared with everyone, so this can take a day.

## 6. While they run (optional)

```bash
./hpc progress gam_lba,gam_lba_ebbinghaus,gam_lba_verticalhorizontal,gam_rdm5_ebbinghaus
```

(B: your four models instead.) This shows, for every running shard, how far
each chain has got and a rough ETA. Once a chain is past warmup (after ~60%
of its time), it also shows the chain's step size and the level of its log
density (`lp__`). There is nothing for you to act on. Pasting the output to
Dom around day 3-4 is enough to spot a broken chain without waiting for the
end.

If a task disappears from `./hpc queue` much earlier than expected, look at
`./hpc log <job number>`. A failed task ends with `status 1`. Resubmit only
that shard, e.g. shard 3 of `gam_lba`:

```bash
./hpc fit gam_lba --array=3 --constraint=amd_zen5
```

Finished shards are never refitted. Tell Dom either way.

## 7. When a model's four shards are done

The model has no task left in `./hpc queue`. For each model:

```bash
./hpc combine gam_lba                              # prints a job number, ~1 h, 320 GB of memory
./hpc predict gam_lba --dependency=afterok:<that job number>
./hpc priors gam_lba --dependency=afterok:<that job number>
```

You can do all four models in one go once the last one is done. When the
jobs have finished (`./hpc queue` is empty), fetch the results:

```bash
./hpc pull                 # combined fits      -> analysis/models/
./hpc pull predictions     # prediction files   -> analysis/models/predictions/
./hpc pull priors          # prior checks       -> analysis/models/priors/
./hpc log <fit job number> > <model>_fitlog.txt    # once per model
```

Then send Dom, through Sussex OneDrive, for each model:

- `analysis/models/<model>_<Illusion>.rds`, the combined fit (~1 GB);
- `analysis/models/predictions/<model>_<Illusion>.rds` (~20 MB);
- `analysis/models/priors/<model>_<Illusion>.rds` (small);
- `<model>_fitlog.txt`.

Do not rename anything: the names are what the notebooks look for. Leave the
shards on the cluster. They are not needed for the transfer, and Dom will say
when they can go.

## If something goes wrong

| symptom | meaning | fix |
| --- | --- | --- |
| `cogmod 0.3.4 is too old; need >= 0.3.5` | step 3 not done | `./hpc install cogmod` |
| `... is from another run (cogmod ...)` | an old shard is still in `models/` | `./hpc archive`, then resubmit |
| compile error on every model, within a minute | the brms setup issue (`AGENT.md` §3.8) | tell Dom |
| `refusing: igc_ jobs are queued or running` | `archive` protects running fits | wait, or `./hpc cancel <job number>` |
| `combine` stops with "different cogmod versions" or "differ in run length" | a shard from another run | tell Dom; nothing is lost |
| `kex_exchange_identification: Connection reset` | too many SSH connections in a row | wait a minute and retry |
