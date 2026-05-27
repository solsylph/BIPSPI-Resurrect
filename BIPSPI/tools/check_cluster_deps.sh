#!/usr/bin/env bash
# Cluster dependency probe for BIPSPI sequence-mode training.
# Run on haskell (or any HPC login node) and pipe stdout back to the chat.
#
# Usage:
#   bash check_cluster_deps.sh
#
# Outputs (sequence-mode requirements, in priority order):
#   1. psiblast        — REQUIRED for PSSM generation. Check `module spider`.
#   2. uniref90 BLAST DB — REQUIRED for psiblast. Check $UNIREF90_DB and ~/Tesis.
#   3. al2co + clustalw + cd-hit — REQUIRED for conservation features.
#   4. SPIDER2 (pred_pssm.py) — REQUIRED for predicted ASA + SS features.
#                              SPIDER2 itself is Python 2.7; flag if found.
#   5. conda/miniforge — required to build the Py3 env.
#   6. xgboost, biopython, sklearn — checked inside the activated env later.
#
# Exit code is 0 always (probe never fails); failures appear as "MISSING" lines.

set -u
echo "=== BIPSPI cluster dependency probe ==="
echo "Date:      $(date -Iseconds)"
echo "Host:      $(hostname)"
echo "User:      $(whoami)"
echo "Shell:     $SHELL"
echo ""

probe_module() {
  local name="$1"
  echo "--- module spider $name ---"
  if command -v module >/dev/null 2>&1; then
    # `module spider` writes to stderr; capture it.
    module spider "$name" 2>&1 | head -30 || echo "MODULE_SPIDER_FAILED"
  else
    echo "MODULE_COMMAND_NOT_AVAILABLE"
  fi
  echo ""
}

probe_binary() {
  local name="$1"
  echo "--- which $name ---"
  if command -v "$name" >/dev/null 2>&1; then
    local path
    path=$(command -v "$name")
    echo "FOUND: $path"
    # Try to get version; many tools take -h, --version, or just exit nonzero.
    { "$path" --version 2>&1 || "$path" -V 2>&1 || true; } | head -5
  else
    echo "MISSING"
  fi
  echo ""
}

probe_file() {
  local desc="$1"
  local path="$2"
  echo "--- $desc: $path ---"
  if [ -e "$path" ]; then
    ls -lah "$path" 2>&1 | head -3
  else
    echo "MISSING"
  fi
  echo ""
}

echo "########## 1. PSI-BLAST ##########"
probe_module psiblast
probe_module BLAST+
probe_module blast
probe_binary psiblast
probe_binary makeblastdb

echo "########## 2. uniref90 BLAST DB ##########"
echo "Looking for a compiled BLAST DB (.phr/.pin/.psq or .pal). Common locations:"
for cand in \
  "$HOME/databases/uniref90.fasta" \
  "$HOME/databases/uniref90" \
  "$HOME/Tesis/databases/NR_SEQ/uniref90/uniref90.blastdb" \
  "$HOME/Tesis/databases/NR_SEQ/uniref90" \
  "/scratch/databases/uniref90" \
  "/data/databases/uniref90" \
  "/shared/databases/uniref90" \
  ; do
  if [ -e "${cand}.pal" ] || [ -e "${cand}.phr" ] || [ -d "$cand" ]; then
    echo "POSSIBLE: $cand"
    ls -lah "$cand"* 2>&1 | head -3
  fi
done
echo ""

echo "########## 3. AL2CO + clustalw + cd-hit ##########"
probe_module al2co
probe_module clustalw
probe_module clustal
probe_module cdhit
probe_module cd-hit
probe_binary al2co
probe_binary clustalw
probe_binary clustalw2
probe_binary cd-hit
probe_binary cdhit

echo "########## 4. SPIDER2 ##########"
echo "SPIDER2 is a Python script bundled with the SPIDER2 distribution."
echo "Searching for pred_pssm.py and SPIDER2 directory:"
for cand in \
  "$HOME/SPIDER2/misc/pred_pssm.py" \
  "$HOME/Tesis/dependencies/bioinformaticTools/SPIDER2/misc/pred_pssm.py" \
  "/shared/tools/SPIDER2/misc/pred_pssm.py" \
  ; do
  probe_file "SPIDER2 script" "$cand"
done
find "$HOME" -maxdepth 4 -name "pred_pssm.py" 2>/dev/null | head -5

echo "########## 5. conda / miniforge ##########"
probe_binary conda
probe_binary mamba
probe_file "miniforge install" "$HOME/miniforge3"
probe_file "active env hint" "$HOME/miniforge3/envs"

echo "########## 6. Python / pip baseline ##########"
probe_binary python3
probe_binary python
echo "--- pip show (in current python) ---"
python3 -c "import sys; print(sys.version)" 2>&1 | head -2
python3 -c "import Bio; print('biopython', Bio.__version__)" 2>&1 | head -2
python3 -c "import xgboost; print('xgboost', xgboost.__version__)" 2>&1 | head -2

echo "########## 7. Resource limits ##########"
echo "--- ulimit -a (head) ---"
ulimit -a 2>&1 | head -10
echo "--- disk free home ---"
df -h "$HOME" 2>&1 | head -3
echo "--- disk free scratch (if any) ---"
df -h /scratch 2>&1 | head -3
df -h /tmp 2>&1 | head -3

echo ""
echo "=== probe complete ==="
echo "Paste this output back to the chat for the implementor to triage."
