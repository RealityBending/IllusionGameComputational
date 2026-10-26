# Onboarding a new cluster account (for the colleague's coding agent)

You are setting up a new lab member's machine and Artemis account (the
University of Sussex HPC) so that they can fit this project's models. The
work ends with a passing smoke test. **Do not submit production fits**: Dom
(Dominique Makowski) decides which models each account runs.

Background reading, in this order, as you need it:

- Lab hub, general to every project: [`hpc/setup.md`](https://github.com/RealityBending/Lab/blob/main/hpc/setup.md)
  (machine, account, project), [`hpc/artemis.md`](https://github.com/RealityBending/Lab/blob/main/hpc/artemis.md)
  (partitions, quotas, `sussexneuro`), [`hpc/toolchain.md`](https://github.com/RealityBending/Lab/blob/main/hpc/toolchain.md)
  (R, CmdStan, the brms trap) and [`hpc/troubleshooting.md`](https://github.com/RealityBending/Lab/blob/main/hpc/troubleshooting.md).
- This project: `analysis/server/README.md` (the `./hpc` driver, and
  "Running from a second cluster account") and `analysis/server/AGENT.md` (the
  project's run history and decisions; long, so search it).

## What only the person can do

Ask them for these, and wait. Do not try to work around them:

1. **An Artemis account.** It is requested from Sussex research computing
   ([Access](https://artemis-docs.hpc.sussex.ac.uk/artemis/access.html)). Note
   their cluster user name (e.g. `abc12`).
2. **The GlobalProtect VPN**, portal `bond.sussex.ac.uk`, with their Sussex
   password and an authenticator code. It must be connected for every step
   that touches the cluster. Without it, DNS resolves but port 22 does not
   answer.
3. **Pasting their public key into the cluster** (step A3 below), through the
   Open OnDemand web shell.
4. **Access to the `sussexneuro` partition**: they ask Dom to have them added
   to the Unix group `artemis_sussexneuro`. It is optional. `long` is enough
   to start.

## A. Their machine (once per laptop)

On Windows, use Git Bash for every command here.

1. Clone both repositories (the project, and the lab hub for its scripts):

   ```bash
   git clone https://github.com/RealityBending/IllusionGameComputational.git
   git clone https://github.com/RealityBending/Lab.git
   ```

2. Create the SSH key and the `artemis` host alias:

   ```bash
   bash Lab/hpc/scripts/setup-ssh.sh <user>
   ```

   It is idempotent. It writes `~/.ssh/artemis` (no passphrase, so you can use
   it without prompts) and an `artemis` entry in `~/.ssh/config`, and prints a
   block to paste into the cluster. **Never copy a private key from another
   machine**: each machine gets its own key.
3. The person opens Open OnDemand → Clusters → `>_ artemis Shell Access` and
   pastes that block. It appends the key to `~/.ssh/authorized_keys`.
4. Check without any chance of a password prompt:

   ```bash
   ssh -o BatchMode=yes -o ConnectTimeout=10 artemis 'whoami; hostname'
   ```

   It should print their user name and `artemis-login-0`. Use this check, not
   `./hpc check`, while the key may not be installed yet: `./hpc` waits on a
   password prompt and hangs instead of failing. If it fails, diagnose with
   `ssh -v -o BatchMode=yes artemis true`.

## B. Their account (once per person)

1. Look at the login environment (hub `setup.md` B1):

   ```bash
   ssh artemis 'cat ~/.bashrc ~/.Renviron 2>/dev/null; id -Gn'
   ```

   They need no `R_LIBS_USER` line: `./hpc` sets it in every job, so that the
   fits get the same brms 2.21.0 as Dom's (`AGENT.md` §3.8). Report any R
   setting you see there anyway. `id -Gn` lists their groups.
   `artemis_sussexneuro` there means they have `sussexneuro`.
2. Point the project driver at their account. `hpc.local` is gitignored, so
   no tracked file changes:

   ```bash
   cd IllusionGameComputational/analysis/server
   echo 'IGC_HPC_USER=<user>' > hpc.local
   ```

   If their lustre directory is not under `psych`
   (`ls /mnt/lustre/users/` on the cluster shows the school names), add
   `IGC_HPC_GROUP=<school>` on a second line.
3. Build everything. `install` and `precompile` run on a compute node, and
   the first `install` builds CmdStan, about 25 minutes. Run it in the
   background:

   ```bash
   ./hpc check        # VPN + SSH
   ./hpc setup        # project directories on lustre
   ./hpc push         # scripts to the cluster (strips CRLF)
   ./hpc install      # R library + CmdStan
   ./hpc precompile   # CmdStan precompiled header
   ```

   `install` must end with every package `OK`, including `brms OK 2.21.0`
   and `cogmod >= 0.3.5 OK`. Note the version on its `cmdstan:` line (Dom's
   account runs 2.39.0).

## C. Smoke test

Fit 30 participants on the 2-hour `short` partition, in a directory of its
own so that it can never mix with real fits:

```bash
IGC_MODELS_DIR=/mnt/lustre/users/psych/<user>/IGComputational/smoke_onboarding \
IGC_NPARTICIPANTS=30 IGC_WARMUP=300 IGC_SAMPLES=100 \
  ./hpc fit gam_lnr_mullerlyer,gam_ddm4_mullerlyer --array=1 --partition=short --cpus-per-task=8 --mem=16G
```

(replace `psych` if step B2 set another group). Then, for each job number
that `sbatch` prints:

```bash
./hpc queue
./hpc log <job number>
```

A pass is `brms 2.21.0 from /research/cisc2/...` near the top of the log,
then `participants: 30  rows: 3734`, a `REPORT` line, `SUCCESSFUL` and
`status 0`, within about 5-15 minutes. Any other brms line, or a compile
error on both models within the first minute, is the brms trap
(`AGENT.md` §3.8): stop and report it.

If they have `sussexneuro`, confirm it with a bare submission. A partition
list is not a test, because `-p sussexneuro,general` silently runs on
`general`:

```bash
./hpc sh "sbatch -p sussexneuro -t 1 -c 1 --mem=100M --wrap=hostname"
```

## Then stop and report

Send Dom (through the person) the user name, the output of `id -Gn`,
whether `sussexneuro` worked, the CmdStan version, and both smoke-test logs. Dom replies with the
models to fit and the exact `./hpc fit` lines. `refit_colleagues.md` shows
what that phase looks like (submit, check progress, combine, send the files
back).

## Rules on this cluster

- **Do not retry a failed SSH connection in a loop.**
  `kex_exchange_identification: Connection reset` means sshd is
  rate-limiting a burst of connections: wait a minute, then try once more.
  A connection that times out means the VPN is down: ask the person.
- Go through `./hpc` for everything. `./hpc sh '<cmd>'` runs a command on the
  login node. Never run R, compiles or fits there: `./hpc` sends them to
  compute nodes.
- Never let a `.slurm` or `.R` file with CRLF line endings reach the cluster.
  `./hpc push` strips them, and `.gitattributes` keeps them out of git.
- Do not reinstall cogmod (`./hpc install cogmod`) while any job of theirs is
  pending or running. A shard that starts after a reinstall is fitted on
  different priors, and the combine step refuses to merge it.
- Ask the person before any `./hpc fit` beyond the smoke test, and before
  `./hpc cancel` or `./hpc archive`.
- The cluster clock is UTC.
