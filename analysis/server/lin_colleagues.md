# Fitting a Study 2 linear model on your own Artemis account (November 2026)

Study 2 fits one model per illusion: `lin_ddm5`, the DDM-5 of Study 1 with
linear effects of difficulty and illusion strength instead of smooth
surfaces (`analysis/3_linearmodels.qmd`). Dom fits the Ebbinghaus one. You are
taking one of the other two:

| model | illusion |
| --- | --- |
| `lin_ddm5_mullerlyer` | Müller-Lyer |
| `lin_ddm5_verticalhorizontal` | Vertical-Horizontal |

Dom's message tells you which is yours (or both: one account can run the two
at once). Below, `<model>` is that name and `<you>` your cluster user name.

What is different from the October refit (`refit_colleagues.md`):

- **The sampler is stanli, not cmdstanr**: stanli 0.31.0, with its own NUTS
  and eight threads per chain, as Dom's Ebbinghaus fit, so that the three
  illusions are sampled the same way. Your account needs its own copy of
  stanli, once (step 2, about 10 minutes): another account's files cannot be
  read.
- **Shorter runs**: 1,000 warmup + 1,000 draws per chain, 8 chains per model
  (4 shards of 2), set on the command line in step 4.
- **cogmod stays where it is.** Do **not** run `./hpc install cogmod`: GitHub's
  `dev` is now 0.3.6, which changes the DDM likelihood, and all three
  illusions must be fitted on the same cogmod, down to the commit (Dom's is
  0.3.5 at `eb4c6a3`). The smoke test shows yours.

All commands run from `analysis/server/` in your clone of the repo, with the
**GlobalProtect VPN connected**. On Windows, use Git Bash. `README.md` in the
same folder is the full reference.

## 1. Update and settings (5 minutes)

```bash
git pull
cd analysis/server
./hpc check            # VPN + SSH
./hpc queue            # see what of yours is running; nothing needs to stop
```

**First time on this project?** Follow README → "Running from a second cluster
account" first. If you fitted models in October, skip it.

Then add three lines to `analysis/server/hpc.local`, the file that already
holds your `IGC_HPC_USER` line (it is gitignored, so nothing you commit
changes):

```bash
IGC_BACKEND=stanli
IGC_STANLI_LIB=/mnt/lustre/users/psych/<you>/cogmod_stanli/R_libs_0.31.0
IGC_STANLI_RUNTIME=/mnt/lustre/users/psych/<you>/cogmod_stanli/runtime_v0.31.0/libstanli.so
```

(`psych` is your school's directory; change it if yours differs, as in
`IGC_HPC_GROUP`.) From then on every `./hpc` command works on stanli fits,
which go to `IGComputational/models_stanli/` on the cluster and to
`analysis/models_stanli/` on your laptop, away from your October fits in
`models/`. **Remove the three lines when you are done** (step 7).

## 2. Install stanli (10-15 minutes, once)

```bash
./hpc push             # copies the current scripts to the cluster
./hpc install-stanli
```

This downloads stanli 0.31.0 and its runtime from the GitHub release, checks
them against their recorded sha256, builds stanli on a compute node of the
`short` partition and tests it. It must end with `STANLI INSTALL OK` and print
the same two paths as your `hpc.local` lines. If `short` is busy, it waits
for a node; `IGC_INSTALL_PARTITION=general ./hpc install-stanli` uses another
partition.

## 3. Smoke test (about 10 minutes)

30 participants, a few iterations, on the 2-hour `short` partition, in a
directory of its own so that it can never be mixed into the real fit:

```bash
IGC_MODELS_DIR=/mnt/lustre/users/psych/<you>/IGComputational/smoke_lin \
IGC_SAMPLER=nuts IGC_CHAINS=2 IGC_WARMUP=100 IGC_SAMPLES=50 \
IGC_NPARTICIPANTS=30 IGC_PARTICIPANT_SEED=2026 \
  ./hpc fit <model> --array=1 --partition=short --cpus-per-task=8 --mem=16G --time=02:00:00
```

Then, with the job number `sbatch` printed:

```bash
./hpc log <job number>
```

You want, in the `.out`:

- `stanli 0.31.0 from /mnt/lustre/users/psych/<you>/cogmod_stanli/R_libs_0.31.0/stanli`;
- **`cogmod 0.3.5 eb4c6a394a72b4e14108fcef00ab5c7c8e107ef3`**, near the top,
  and the same at the end (`cogmod of this shard: ...`). Any other version
  or commit: stop and tell Dom before step 4;
- `participants: 30`;
- `stanli inits: one-element arrays: sd_...`, the fix for what stopped Dom's
  first Ebbinghaus run;
- a `REPORT` line and `FINISHED ... with status 0`. On Dom's account it took
  9 minutes. Its step size, tree depth and Rhat are poor after only 100
  warmup iterations, which is expected and says nothing about the real fit.

## 4. Submit the real fit

```bash
IGC_SAMPLER=nuts IGC_CHAINS=2 IGC_WARMUP=1000 IGC_SAMPLES=1000 \
  ./hpc fit <model> --array=1-4 --partition=long --constraint=amd_zen5 \
    --cpus-per-task=16 --mem=64G --time=8-00:00:00
```

(Both models: `./hpc fit lin_ddm5_mullerlyer,lin_ddm5_verticalhorizontal`
with the same options. Eight tasks of 16 CPUs fit in `long`'s 140-CPU cap.)

Use exactly that: `long` (not `sussexneuro`), `--constraint=amd_zen5`, 64G,
8 days. **Note the job number.** `./hpc queue` should show four tasks,
probably `PENDING (Resources)` or `(Priority)` at first: they wait for a zen5
node, of which there are three, and that can take a day.

How long a shard takes at full data has not been measured for this model.
The estimate is half a day to two days once it has a node. The 8 days are
counted per shard from its start, and are a margin. A stanli shard saves
nothing until it ends, so a shard that hits the limit is lost.

## 5. While it runs (optional)

```bash
./hpc progress <model>
```

This shows each chain's iteration and a rough ETA. stanli writes no CmdStan
CSV, so unlike the October fits there is no step size or `lp__` here.
Pasting it to Dom once a day is enough. A task that disappears from
`./hpc queue` early has failed: `./hpc log <job number>`, then tell Dom.
To rerun one shard, e.g. shard 3, repeat step 4 with `--array=3`: finished
shards are never refitted.

## 6. When the four shards are done

```bash
./hpc combine <model>      # prints a job number; adds loo, about an hour
./hpc predict <model> --dependency=afterok:<that job number> --partition=short
```

When both have finished (`./hpc queue` is empty):

```bash
./hpc pull predictions/<model>_<Illusion>.rds     # -> analysis/models_stanli/predictions/
./hpc log <fit job number> > <model>_fitlog.txt
```

with `<Illusion>` `MullerLyer` or `VerticalHorizontal`, e.g.
`predictions/lin_ddm5_mullerlyer_MullerLyer.rds`. Send Dom, through Sussex
OneDrive, the prediction file (tens of MB) and the fit log. The analyses read
only the prediction file. Leave the combined fit and the shards on the
cluster: do not pull them, and Dom will say when they can go.

## 7. Afterwards

Remove the three lines of step 1 from `hpc.local`. Keep `cogmod_stanli/` and
the fits until Dom says otherwise.

## If something goes wrong

| symptom | meaning | fix |
| --- | --- | --- |
| `refusing: stanli 0.31.0 is already installed` | step 2 was done before | nothing to do |
| `sha256sum: WARNING: ... did NOT match` | the download is not the release | tell Dom; do not build it |
| `there is no package called 'stanli'`, or a runtime error at load | the `hpc.local` paths do not match the install | compare them with what step 2 printed |
| `mismatch in number dimensions ... dims found=()` | your scripts predate the fix | `git pull`, `./hpc push`, resubmit |
| another cogmod version or commit than `0.3.5 eb4c6a3...` in the log | your library was installed at another commit of `dev` | tell Dom before submitting. The fix, once Dom agrees and nothing of yours is queued, is `IGC_COGMOD_REF=eb4c6a394a72b4e14108fcef00ab5c7c8e107ef3 ./hpc install cogmod` |
| `cogmod ... is too old; need >= 0.3.5` | your library predates the October refit | tell Dom |
| `kex_exchange_identification: Connection reset` | too many SSH connections in a row | wait a minute and retry |
