# Node history: which Artemis nodes were fast and which were slow

A running log of every full-data fit shard, the node it landed on, and how long
it took, to find out whether node placement predicts wall time. AGENT.md 4.4
measured a 2-3x per-gradient difference between an idle zen5 node and a loaded
zen3 one on a 120-participant fit. This file tracks the same question at
production scale. **Add a row whenever a shard finishes or is cancelled.**

All rows: full data (2,215 participants, MullerLyer), warmup 1000 +
500 samples, 2 chains x 8 threads per shard, cold start. "Wall" is the
`REPORT` line's wall time for the shard (both chains run in parallel). For a
shard that did not finish, it is how far it got and after how long.

Node CPUs (from `scontrol show node`):
- `artemis-rtx-00/01/02`: **amd_zen5**, 128 CPUs
- `artemis-a40-02/11`: **amd_zen3**, 128 CPUs (a40-11 is also highmem)
- `artemis-cpu-01`: **amd_zen4**, CPU-only, highmem
- `artemis-a40-16`: **amd_zen3**, highmem
- `artemis-general-02`: `virtual`, 128 CPUs (CPU type not advertised)

The rtx nodes are the only zen5 nodes in `sussexneuro` (`sinfo -N -p sussexneuro -o "%N %f"`), so `--constraint=amd_zen5` means a pool of three.

## Log

| model | job | shard | node | CPU | outcome | wall |
|---|---|---|---|---|---|---|
| gam_lnr | 11385713 | 1 | a40-16 | zen3 | cancelled at 600/1500 | 47.0 h |
| gam_lnr | 11385713 | 2 | cpu-01 | zen4 | done | 44.0 h |
| gam_lnr | 11385713 | 3 | cpu-01 | zen4 | done | 42.7 h |
| gam_lnr | 11385713 | 4 | cpu-01 | zen4 | done | 42.8 h |
| gam_lnr6 | 11385717 | 1 | a40-16 | zen3 | cancelled at 300/1500 | 46.9 h |
| gam_lnr6 | 11385717 | 2-4 | cpu-01 | zen4 | cancelled at 600-800/1500 | 46.9 h |
| gam_lnr | 11399599 | 1 | rtx-02 | zen5 | done | 25.5 h |
| gam_lnr | 11399599 | 2 | rtx-02 | zen5 | done | 29.0 h |
| gam_lnr | 11399599 | 3 | rtx-02 | zen5 | done | 16.9 h |
| gam_lnr | 11399599 | 4 | rtx-01 | zen5 | done | 22.5 h |
| gam_lnr6 | 11399603 | 1 | rtx-01 | zen5 | done | 72.8 h |
| gam_lnr6 | 11399603 | 2 | rtx-01 | zen5 | done | 70.9 h |
| gam_lnr6 | 11399603 | 3 | rtx-01 | zen5 | done | 85.9 h |
| gam_lnr6 | 11399603 | 4 | a40-11 | zen3 | **cancelled at 500/1500** | 108 h |
| gam_ddm4 | 11403693 | 2 | rtx-00 | zen5 | done | 25.3 h |
| gam_ddm4 | 11403693 | 3 | rtx-00 | zen5 | done | 38.6 h |
| gam_ddm4 | 11403693 | 4 | rtx-00 | zen5 | done | 32.2 h |
| gam_ddm4 | 11403693 | 1 | general-02 | virtual | still running at 1300/1500 | 72.7 h |
| gam_ddm5 | 11404322 | 1 | rtx-02 | zen5 | done | 33.5 h |
| gam_ddm5 | 11404322 | 2-4 | a40-11 | zen3 | still running at 700-1200/1500 | 59.3 h |
| gam_lnr6 | 11406608 | 5 | a40-02 | zen3 | **cancelled at 200/1500** | 38.6 h |
| gam_lnr6 | 11410717 | 6 | rtx-01 | zen5 (`--constraint`) | running, started 2026-09-24 | |

## What it says so far (2026-09-24)

- **Every zen5 (`rtx`) shard finished, and every slow shard was on something
  else.** The cleanest evidence compares shards *within the same array*, which
  share the code, data and submission time and differ only in the node:
  - gam_lnr6 11399603: 71-86 h on rtx-01; the a40-11 shard had done a third
    of its iterations after 108 h, so it was roughly 4x slower.
  - gam_ddm5 11404322: 33.5 h on rtx-02; after 59 h the a40-11 shards are
    still 47-80% of the way through.
  - gam_ddm4 11403693: 25-39 h on rtx-00; general-02 took about 2x as long.
- The backup gam_lnr6 shard on a40-02 (zen3, lightly loaded) was as slow as
  the one on a40-11, so the slowdown looks like the CPU generation rather than
  one overloaded node.
- The 2026-09-18 runs (cpu-01 and a40-16) predate the 2026-09-20 rerun, so
  comparing them with the rtx rows mixes two code versions. Within that run,
  cpu-01 (zen4) finished gam_lnr in about 43 h, while a40-16 (zen3) had done
  40% after 47 h. Across all runs the order is zen5 > zen4 > zen3.
- **Caveat:** so far that is a handful of shards per node type, and node load
  varies (rtx-01 carried 94 of its 128 CPUs and was still fast). This is a
  strong pattern, not a controlled measurement.

## Rule of thumb until the evidence says otherwise

Pin long fits to zen5 nodes:

```bash
./hpc fit <model> --constraint=amd_zen5
```

The cost is a smaller pool of nodes, so a job can queue for longer. Check
`squeue` after submitting. If it sits in `PENDING` with reason `Resources`
while zen3 nodes are free, the wait may cost more than the slowdown would.
