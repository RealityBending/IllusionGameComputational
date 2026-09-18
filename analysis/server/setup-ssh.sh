#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# One-shot Artemis SSH setup for a NEW machine.
#
#   bash analysis/server/setup-ssh.sh
#
# Generates a machine-local keypair (never copy private keys between
# machines), adds the 'artemis' host alias to ~/.ssh/config, and prints the
# public key plus the exact command to install it on the cluster.
# Idempotent: safe to re-run.
# ---------------------------------------------------------------------------
set -euo pipefail

HPC_USER="${IGC_HPC_USER:-dmm56}"
HPC_HOST=ood.artemis.hrc.sussex.ac.uk
KEY="$HOME/.ssh/artemis"

mkdir -p "$HOME/.ssh" && chmod 700 "$HOME/.ssh"

if [ -f "$KEY" ]; then
  echo "==> key already exists: $KEY"
else
  echo "==> generating ed25519 keypair for $(hostname)"
  ssh-keygen -t ed25519 -f "$KEY" -N "" -C "${HPC_USER}@$(hostname)"
  ssh-keygen -t rsa -b 4096 -f "${KEY}_rsa" -N "" -C "${HPC_USER}@$(hostname)-rsa"
fi
chmod 600 "$KEY" "${KEY}_rsa" 2>/dev/null || true

CFG="$HOME/.ssh/config"
if [ -f "$CFG" ] && grep -q "^Host artemis$" "$CFG"; then
  echo "==> ~/.ssh/config already has the 'artemis' alias"
else
  echo "==> adding 'artemis' alias to ~/.ssh/config"
  cat >> "$CFG" <<EOF

# --- Artemis HPC (University of Sussex) ---------------------------------
# Requires the GlobalProtect VPN (portal bond.sussex.ac.uk) to be connected.
Host artemis
    HostName ${HPC_HOST}
    User ${HPC_USER}
    IdentityFile ~/.ssh/artemis
    IdentitiesOnly yes
    StrictHostKeyChecking accept-new
    ConnectTimeout 15
    ServerAliveInterval 60
    ServerAliveCountMax 3

Host artemis-rsa
    HostName ${HPC_HOST}
    User ${HPC_USER}
    IdentityFile ~/.ssh/artemis_rsa
    IdentitiesOnly yes
    StrictHostKeyChecking accept-new
    ConnectTimeout 15
EOF
  chmod 600 "$CFG"
fi

cat <<EOF

---------------------------------------------------------------------------
Next: install this machine's public key on Artemis.

With the VPN connected, open
  https://ood.artemis.hrc.sussex.ac.uk/
then Clusters -> ">_ artemis Shell Access" and paste this whole block:

mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys <<'KEYS'
$(cat "${KEY}.pub")
KEYS
chmod 600 ~/.ssh/authorized_keys && echo OK

Then check it from here with:
  cd analysis/server && ./hpc check
---------------------------------------------------------------------------
EOF
