#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_remaining.sh -- send the fits that were still missing on 2026-09-28:
# the MullerLyer models on VerticalHorizontal and Ebbinghaus (all defined in
# models.R). Run the stages in order, from analysis/server, with the VPN up,
# after committing models.R and running `./hpc push` (the jobs read the
# registry on the cluster, not this copy):
#
#   ./submit_remaining.sh smoke     1. Ebbinghaus has never been fitted: one
#                                      30-participant LNR on `short`, in its own
#                                      directory. Check `./hpc log` before 2.
#   ./submit_remaining.sh long      2. the production queue on `long`, zen5 only
#   ./submit_remaining.sh neuro     3. the slowest model on `sussexneuro`, which
#                                      does not draw on `long`'s quota
#   ./submit_remaining.sh lba       4. the two LBAs -- only once the MullerLyer
#                                      gam_lba (a colleague's run) has a clean
#                                      REPORT line
#
# Add --test-only to any stage to have SLURM validate the requests and say
# where they would start, without queueing anything.
#
# Each stage is one `./hpc` call, i.e. two SSH connections however many models
# it submits (./hpc fit takes a comma-separated list).
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
[ -f hpc.local ] && . ./hpc.local

stage="${1:-}"; shift || true

# Submission order is run order: on `long` one account runs 8 tasks (two
# models) at a time and the rest wait, in the order submitted, at no cost.
# VerticalHorizontal first, since it only lacks these four; the DDM-5 is left
# out here because it goes to sussexneuro (stage 3).
LONG=(
  gam_lnr6_verticalhorizontal
  gam_rdm_verticalhorizontal
  gam_rdm5_verticalhorizontal
  gam_lnr_ebbinghaus
  gam_ddm4_ebbinghaus
  gam_rdm_ebbinghaus
  gam_lnr6_ebbinghaus
  gam_rdm5_ebbinghaus
  gam_ddm5_ebbinghaus
)
# The slowest model (DDM-5: 33.5 h per shard on zen5, against 15-27 h for the
# others), so it gets the partition that lets it start now.
NEURO=gam_ddm5_verticalhorizontal
# Held back: gam_lba has not finished on MullerLyer yet, and at warmup 300 its
# smoke test had 6% divergences (AGENT.md 4.7). Fitting it twice more before
# seeing its first production REPORT would risk 8 of the queue's shards.
LBA=(gam_lba_verticalhorizontal gam_lba_ebbinghaus)

join() { local IFS=,; echo "$*"; }
USERS_DIR="${IGC_USERS_DIR:-/mnt/lustre/users/${IGC_HPC_GROUP:-psych}/${IGC_HPC_USER:-dmm56}/${IGC_PROJECT:-IGComputational}}"

case "$stage" in
  smoke)
    # Its own IGC_MODELS_DIR: fits use file_refit = "never", so a 30-person
    # shard left in the production directory would be adopted by stage 2.
    IGC_MODELS_DIR="${USERS_DIR}/smoke_ebbinghaus" \
    IGC_NPARTICIPANTS=30 IGC_WARMUP=300 IGC_SAMPLES=100 \
      ./hpc fit gam_lnr_ebbinghaus --array=1 --partition=short --cpus-per-task=8 --mem=16G "$@"
    echo "Then: ./hpc log gam_lnr_ebbinghaus -- want a 'participants: 30 rows: ~4300' line,"
    echo "a REPORT line and 'status 0' (the MullerLyer LNR took 6-11 min)."
    ;;
  long)
    ./hpc fit "$(join "${LONG[@]}")" --constraint=amd_zen5 "$@"
    ;;
  neuro)
    # sussexneuro's 256 CPUs are one pool for the whole group, and a request
    # past it is rejected rather than queued: look before asking for 64.
    used=$(./hpc sh "squeue -p sussexneuro -h -t RUNNING,PENDING -o %C" 2>/dev/null |
           awk '/^[0-9]+$/ {s += $1} END {print s + 0}')
    echo "sussexneuro: ${used} of 256 group CPUs in use or requested"
    if [ "$used" -gt $((256 - 64)) ]; then
      echo "not enough room for 4 x 16 CPUs -- run 'long' with ${NEURO} added instead,"
      echo "or retry later." >&2
      exit 1
    fi
    ./hpc fit "$NEURO" --partition=sussexneuro --constraint=amd_zen5 "$@"
    ;;
  lba)
    ./hpc fit "$(join "${LBA[@]}")" --constraint=amd_zen5 "$@"
    ;;
  *)
    sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 1
    ;;
esac
