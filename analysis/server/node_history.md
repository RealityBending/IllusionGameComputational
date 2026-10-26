# Node history: which Artemis nodes were fast and which were slow

A running log of every full-data fit shard, the node it landed on, and how long
it took, to find out whether node placement predicts wall time. AGENT.md 4.4
measured a 2-3x per-gradient difference between an idle zen5 node and a loaded
zen3 one on a 120-participant fit. This file tracks the same question at
production scale. **Add a row whenever a shard finishes or is cancelled.**

All rows up to 2026-10-13: full data, warmup 1000 + 500 samples, 2 chains x 8
threads per shard, cold start; rows from the 2026-10-13 refit onwards are
1500 + 1000 (AGENT.md 4.11) and say so. Each `.out` also records every chain's
warmup and sampling seconds (`REPORT time per chain`), which is where the
per-model costs in AGENT.md 4.11 come from. "Wall" is the
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
| gam_ddm4 | 11403693 | 1 | general-02 | virtual | done | 82.6 h |
| gam_ddm5 | 11404322 | 1 | rtx-02 | zen5 | done | 33.5 h |
| gam_ddm5 | 11404322 | 3 | a40-11 | zen3 | done | 76.6 h |
| gam_ddm5 | 11404322 | 2 | a40-11 | zen3 | done | 116.5 h |
| gam_ddm5 | 11404322 | 4 | a40-11 | zen3 | done | 137.1 h |
| gam_lnr6 | 11406608 | 5 | a40-02 | zen3 | **cancelled at 200/1500** | 38.6 h |
| gam_lnr6 | 11410717 | 6 | rtx-01 | zen5 (`--constraint`) | done, **chain 1 stuck** (step size 2.9e-05, lp__ -55,700 vs +38,900) | 84.8 h |
| gam_lnr_verticalhorizontal | 11411625 | 1 | rtx-00 | zen5 (`--constraint`) | done | 17.5 h |
| gam_lnr_verticalhorizontal | 11411625 | 2 | rtx-00 | zen5 (`--constraint`) | done | 20.6 h |
| gam_lnr_verticalhorizontal | 11411625 | 3 | rtx-00 | zen5 (`--constraint`) | done | 20.7 h |
| gam_lnr_verticalhorizontal | 11411625 | 4 | rtx-02 | zen5 (`--constraint`) | done | 14.5 h |
| gam_ddm4_verticalhorizontal | 11414141 | 1 | rtx-00 | zen5 (`--constraint`) | done | 26.7 h |
| gam_ddm4_verticalhorizontal | 11414141 | 2 | rtx-01 | zen5 (`--constraint`) | done | 17.8 h |
| gam_ddm4_verticalhorizontal | 11414141 | 3 | rtx-01 | zen5 (`--constraint`) | done | 21.4 h |
| gam_ddm4_verticalhorizontal | 11414141 | 4 | rtx-01 | zen5 (`--constraint`) | done | 20.1 h |
| gam_lnr6_verticalhorizontal | 11415945 | 1 | rtx-00 | zen5 (`--constraint`) | done | 41.5 h |
| gam_lnr6_verticalhorizontal | 11415945 | 2 | rtx-01 | zen5 (`--constraint`) | done | 54.2 h |
| gam_lnr6_verticalhorizontal | 11415945 | 3 | rtx-01 | zen5 (`--constraint`) | done | 48.6 h |
| gam_lnr6_verticalhorizontal | 11415945 | 4 | rtx-02 | zen5 (`--constraint`) | done | 57.0 h |
| gam_ddm5_verticalhorizontal | 11415949 | 1 | rtx-02 | zen5 (`--constraint` on sussexneuro) | done | 23.9 h |
| gam_ddm5_verticalhorizontal | 11415949 | 2 | rtx-02 | zen5 (`--constraint` on sussexneuro) | done | 23.8 h |
| gam_ddm5_verticalhorizontal | 11415949 | 3 | rtx-02 | zen5 (`--constraint` on sussexneuro) | done | 22.1 h |
| gam_ddm5_verticalhorizontal | 11415949 | 4 | rtx-02 | zen5 (`--constraint` on sussexneuro) | done | 23.4 h |
| gam_rdm_verticalhorizontal | 11417276 | 1 | rtx-00 | zen5 (`--constraint`) | done | 31.3 h |
| gam_rdm_verticalhorizontal | 11417276 | 2 | rtx-00 | zen5 (`--constraint`) | done | 29.9 h |
| gam_rdm_verticalhorizontal | 11417276 | 3 | rtx-00 | zen5 (`--constraint`) | done | 30.1 h |
| gam_rdm_verticalhorizontal | 11417276 | 4 | rtx-02 | zen5 (`--constraint`) | done | 35.5 h |
| gam_lnr6 | 11417300 | 6 (refit, `IGC_FILE_REFIT=always`) | rtx-02 | zen5 (`--constraint` on sussexneuro) | done, both chains healthy (step 0.005-0.008) | 76.7 h |
| gam_lnr_ebbinghaus | 11419150 | 1 | rtx-02 | zen5 (`--constraint`) | done, **adaptation failed** (step 1e-8 to 3e-5, Rhat 3.05; cogmod 0.3.3) | 35.9 h |
| gam_lnr_ebbinghaus | 11419150 | 2-4 | rtx-01 | zen5 (`--constraint`) | done, **adaptation failed** (step 5e-5 to 1e-3) | 32.3 / 36.5 / 32.0 h |
| gam_rdm5_verticalhorizontal | 11419122 | 1-4 | rtx-00 | zen5 (`--constraint`) | done | 65.6 / 66.7 / 64.6 / 41.1 h |
| gam_ddm4_ebbinghaus | 11423962 | 1 | rtx-00 | zen5 | done | 20.8 h (sacct) |
| gam_ddm4_ebbinghaus | 11423962 | 2 | rtx-01 | zen5 | done | 18.7 h (sacct) |
| gam_ddm4_ebbinghaus | 11423962 | 3 | rtx-00 | zen5 | done | 19.5 h (sacct) |
| gam_ddm4_ebbinghaus | 11423962 | 4 | rtx-01 | zen5 | done | 36.5 h (sacct) |
| gam_rdm_ebbinghaus | 11423963 | 1 | rtx-02 | zen5 | sampled, **saveRDS failed** (truncated shard deleted 2026-10-06) | 23.5 h (sacct) |
| gam_rdm_ebbinghaus | 11423963 | 2 | rtx-02 | zen5 | done | 33.7 h (sacct) |
| gam_rdm_ebbinghaus | 11423963 | 3 | rtx-02 | zen5 | done | 22.3 h (sacct) |
| gam_rdm_ebbinghaus | 11423963 | 4 | rtx-02 | zen5 | sampled, **saveRDS failed** (truncated shard deleted 2026-10-06) | 24.2 h (sacct) |
| gam_ddm5_ebbinghaus | 11423966 | 1 | rtx-01 | zen5 (sussexneuro) | sampled, **saveRDS failed** (truncated shard deleted 2026-10-06) | 32.4 h (sacct) |
| gam_ddm5_ebbinghaus | 11423966 | 2 | rtx-01 | zen5 (sussexneuro) | sampled, **saveRDS failed** (truncated shard deleted 2026-10-06) | 33.3 h (sacct) |
| gam_ddm5_ebbinghaus | 11423966 | 3 | rtx-01 | zen5 (sussexneuro) | sampled, **saveRDS failed** (truncated shard deleted 2026-10-06) | 32.4 h (sacct) |
| gam_ddm5_ebbinghaus | 11423966 | 4 | rtx-01 | zen5 (sussexneuro) | done | 33.6 h (sacct) |
| gam_rdm_ebbinghaus | 11430274 | 1 (refit of a failed save) | rtx-00 | zen5 (`--constraint`) | done | 25.3 h (sacct) |
| gam_rdm_ebbinghaus | 11430274 | 4 (refit of a failed save) | rtx-02 | zen5 (`--constraint`) | done | 30.3 h (sacct) |
| gam_ddm5_ebbinghaus | 11430275 | 1-3 (refit of the failed saves) | rtx-02 | zen5 (`--constraint` on sussexneuro) | done | 27.9 / 24.2 / 26.4 h (sacct) |
| gam_rdm5_ebbinghaus | 11553382 | 1-3 | rtx-01 | zen5 (`--constraint`) | done, converged (AGENT.md 4.10) | 47.4 / 59.6 / 64.0 h |
| gam_rdm5_ebbinghaus | 11553382 | 4 | rtx-00 | zen5 (`--constraint`) | done | 79.4 h |
| gam_lnr6_ebbinghaus | 11553384 | 2, 3 | rtx-02 | zen5 (`--constraint`) | done, **adaptation failed** (step 1e-4 to 5e-4, Rhat 3.05) | 84.5 / 82.0 h |
| gam_lnr6_ebbinghaus | 11553384 | 1, 4 | rtx-02 | zen5 (`--constraint`) | done 2026-10-13, same failure (step 8e-4 to 9e-6, Rhat 3.04-3.05) | 92.8 / 92.8 h (sacct) |
| gam_lnr_ebbinghaus | 11553388 | 1-4 (refit on cogmod 0.3.4) | rtx-00 | zen5 (`--constraint` on sussexneuro) | done, **adaptation failed** in all 8 chains (step 1.4e-4 to 5.7e-4, Rhat 3.04) | 61.9 / 40.7 / 35.0 / 41.6 h |
| gam_ddm4_verticalhorizontal | 11596775 | 1-4 (2026-10-13 refit: 1500 + 1000, cogmod 0.3.5) | rtx-01 | zen5 (`--constraint`) | done 2026-10-18: step 0.016-0.019 in all 8 chains, n_leapfrog 254, divergent 0-0.5%, max Rhat 1.83 per shard. Per chain ~25-28 h warmup + ~14-15 h sampling | 42.3 / 42.7 / 41.9 / 43.2 h |
| gam_ddm4_ebbinghaus | 11596774 | 4 (2026-10-13 refit: 1500 + 1000, cogmod 0.3.5) | rtx-01 | zen5 (`--constraint`) | done 2026-10-18: step 0.0163 / 0.0156, n_leapfrog 255, divergent 0, max Rhat 1.83. Chains 41.4 / 43.1 h | 43.2 h |
| gam_ddm4_ebbinghaus | 11596774 | 2, 3 (same run) | rtx-02 | zen5 (`--constraint`) | done 2026-10-19: step 0.014-0.017, n_leapfrog 254-255, divergent 0.2-0.4%, max Rhat 1.83. Chains 66.1 / 67.1 h and 58.2 / 60.5 h | 67.2 / 60.5 h |
| gam_ddm4_ebbinghaus | 11596774 | 1 (same run) | rtx-00 | zen5 (`--constraint`) | done 2026-10-19 12:51: step 0.0124 / 0.0175, mean n_leapfrog 382 (one chain at tree depth 9), divergent 0.4%, max Rhat 1.95. Chains 56.5 / 72.2 h | 72.3 h |
| gam_lnr_mullerlyer, **stanli** | 11603743 | 1 (stanli speed test, AGENT.md 6.6; 1500 + 1000, cogmod 0.3.5) | rtx-02 | zen5 (`--constraint` on sussexneuro) | done 2026-10-20 05:18: step 0.0064, n_leapfrog 766, treedepth 9.50, divergent 0.3%, max Rhat 1.83. Chains 64.4 / 57.4 h (warmup 45.0 / 44.5 h). Its 4 cmdstanr counterparts (11603736), on the same node, were at 1500-2200 / 2500 at that point | 64.8 h |
| gam_lnr_mullerlyer, **stanli** | 11603743 | 4 (same run) | rtx-00 | zen5 (`--constraint` on sussexneuro) | done 2026-10-20 04:29: step 0.00614, n_leapfrog 767, treedepth 9.50, divergent 0, max Rhat 1.83. Chains 63.6 / 55.3 h (warmup 37.9 / 39.2 h) | 64.0 h |
| gam_lnr_mullerlyer, **stanli** | 11603743 | 3 (same run) | rtx-00 | zen5 (`--constraint` on sussexneuro) | done 2026-10-20 22:46: step 0.00476, n_leapfrog 1023, **treedepth 10.00** in both chains, divergent 0, max Rhat 1.56. Chains 81.9 / 80.3 h (warmup 39.9 / 38.0 h, sampling 42 h each) | 82.3 h |
| gam_lnr_mullerlyer | 11603736 | 2, 4 (cmdstanr arm of the stanli test, 1500 + 1000, cogmod 0.3.5) | rtx-02 | zen5 (`--constraint` on sussexneuro) | done 2026-10-21 01:20 / 01:49: step 0.00683 / 0.00695, n_leapfrog 511 / 632, treedepth 9.00 / 9.13, divergent 0.1 / 0.2%, max Rhat 1.83 | 84.9 / 85.3 h |
| gam_lnr_mullerlyer | 11603736 | 3 (same run) | rtx-02 | zen5 (`--constraint` on sussexneuro) | done 2026-10-21 09:30: step 0.00642, n_leapfrog 766, treedepth 9.49, divergent 0.4%, max Rhat 1.83. Shard 1 was at 2000 / 2500 at 102 h, at treedepth 10 | 93.0 h |
| gam_lnr_mullerlyer | 11603736 | 1 (same run) | rtx-02 | zen5 (`--constraint` on sussexneuro) | done 2026-10-23 00:45: both chains at max treedepth on 99-100% of draws (AGENT.md 6.5). REPORT line not read | 132.3 h (sacct) |
| gam_lba_verticalhorizontal | 11596852 (`oc236`) | 1, 2 (2026-10-13 refit: 1500 + 1000, cogmod 0.3.5) | rtx-01 | zen5 (`--constraint`) | done 2026-10-19 18:10 / 22:07 | 77.6 / 81.6 h (sacct) |
| gam_lba_verticalhorizontal | 11596852 (`oc236`) | 3, 4 (same run) | rtx-01 / rtx-00 | zen5 (`--constraint`) | done 2026-10-20 17:30 / 2026-10-21 14:47 | 101.0 / 122.2 h (sacct) |
| gam_rdm5_ebbinghaus | 11596853 (`oc236`) | 2, 3, 4 (2026-10-13 refit: 1500 + 1000, cogmod 0.3.5) | rtx-00 | zen5 (`--constraint`) | done 2026-10-21 15:27 / 08:29 / 07:11. Shard 1 (rtx-00) still running at 125.9 h | 122.9 / 115.9 / 114.6 h (sacct) |
| gam_lnr6_mullerlyer | 11601935 (`asf25`) | 3 (2026-10-13 refit: 1500 + 1000, cogmod 0.3.5) | rtx-02 | zen5 (`--constraint`) | done 2026-10-20 19:51. Shards 1, 2, 4 (rtx-02) still running at 125.6 h | 103.0 h (sacct) |

## What it says so far (updated 2026-09-26)

- **Every zen5 (`rtx`) shard finished, and every slow shard was on something
  else.** The cleanest evidence compares shards *within the same array*, which
  share the code, data and submission time and differ only in the node:
  - gam_lnr6 11399603: 71-86 h on rtx-01; the a40-11 shard had done a third
    of its iterations after 108 h, so it was roughly 4x slower.
  - gam_ddm5 11404322: 33.5 h on rtx-02 against 76.6 h for the first a40-11
    shard to finish (2.3x); the other two a40-11 shards were
    two-thirds of the way through after 80 h.
  - gam_ddm4 11403693: 25-39 h on rtx-00 against 82.6 h on general-02 (2-3x).
  - gam_lnr6 shard 6 (rtx-01, zen5) reached iteration 300 in 21 h and 600
    in 34 h. Shard 4 (a40-11, zen3) took about 56 h to reach 300, and shard 5
    (a40-02) 17 h for its first 100.
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
