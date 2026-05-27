#!/usr/bin/env bash
# Phase B: install al2co + clustalw.
#
# clustalw is installed via bioconda (version 2.1).  BIPSPI docs spec'd
# clustalw 1.83 but al2co only uses standard MSA commands that 1.83 and
# 2.1 share -- backward compatible.
#
# al2co is a single C file, built from the TheApacheCats/al2co GitHub
# mirror (the original ftp://iole.swmed.edu source is dead).  Per BIPSPI
# docs/repo_help.md we patch the name-buffer size from char[500] to
# char[1024] before compiling.
#
# Outputs:
#   - clustalw binary on $PATH inside the activated env (e.g. ~/miniforge3/envs/protein/bin/clustalw)
#   - ~/tools/al2co/al2co  binary
#
# After this script: edit configFiles/cmdTool/dependencies.cfg to point
# clustalW_path and al2coBin_path at the locations printed at the end.

set -euo pipefail

TOOLS=$HOME/tools
mkdir -p "$TOOLS"

# ---------------------------------------------------------------------------
# 1. clustalw via bioconda
# ---------------------------------------------------------------------------
if ! command -v clustalw >/dev/null 2>&1; then
  echo "=== Installing clustalw via bioconda ==="
  conda install -n protein -c bioconda -c conda-forge clustalw -y
else
  echo "clustalw already on PATH: $(command -v clustalw) -- skipping"
fi

CLUSTALW_PATH=$(conda run -n protein which clustalw 2>/dev/null || command -v clustalw)
if [ -z "$CLUSTALW_PATH" ] || [ ! -x "$CLUSTALW_PATH" ]; then
  echo "ERROR: clustalw install reported success but binary not found"
  exit 1
fi
echo "clustalw: $CLUSTALW_PATH"

# ---------------------------------------------------------------------------
# 2. al2co from GitHub mirror, with buffer patch
# ---------------------------------------------------------------------------
AL2CO_DIR=$TOOLS/al2co
if [ ! -x "$AL2CO_DIR/al2co" ]; then
  echo "=== Cloning al2co from TheApacheCats/al2co ==="
  rm -rf "$AL2CO_DIR"
  git clone https://github.com/TheApacheCats/al2co.git "$AL2CO_DIR"
  cd "$AL2CO_DIR"

  # Patch: enlarge name buffers from 500 to 1024 (per BIPSPI docs/repo_help.md).
  echo "Patching name-buffer sizes 500 -> 1024..."
  shopt -s nullglob
  for f in *.c *.h; do
    [ -f "$f" ] || continue
    sed -i -E 's/\bchar([[:space:]]+[A-Za-z_][A-Za-z0-9_]*)\[500\]/char\1[1024]/g' "$f"
  done
  shopt -u nullglob

  # Build: single C file, standard libm link.
  echo "Building al2co..."
  gcc al2co.c -o al2co -lm

  if [ ! -x ./al2co ]; then
    echo "ERROR: al2co build failed -- binary not produced"
    exit 1
  fi
  cd "$TOOLS"
else
  echo "al2co already present at $AL2CO_DIR/al2co -- skipping"
fi

# ---------------------------------------------------------------------------
# 3. Smoke test al2co binary (no args -> usage)
# ---------------------------------------------------------------------------
echo ""
echo "=== Smoke tests ==="
echo "--- clustalw -version ---"
"$CLUSTALW_PATH" -version 2>&1 | head -5 || true
echo ""
echo "--- al2co (no args, expect usage) ---"
"$AL2CO_DIR/al2co" 2>&1 | head -10 || true
echo ""

# ---------------------------------------------------------------------------
# 4. Print the dependencies.cfg lines you need
# ---------------------------------------------------------------------------
cat <<EOF

=== DONE ===

Update ~/BIPSPI-Resurrect/BIPSPI/configFiles/cmdTool/dependencies.cfg so
these lines point at the binaries (replace the existing values):

    al2coBin_path     $AL2CO_DIR/al2co
    clustalW_path     $CLUSTALW_PATH
    cdHitBin_path     $(command -v cd-hit 2>/dev/null || echo "cd-hit  # after: module load CD-HIT/4.8.1-GCC-12.2.0")
    psiBlastBin       $(command -v psiblast 2>/dev/null || echo "psiblast  # after: module load BLAST+/2.14.1-gompi-2023a")
    psiBlastDB_path   /home/biostruct01/databases/uniref90/uniref90.fasta   # after Phase D

(cd-hit and psiblast come from \`module load\`; keep them as bare names
or use absolute paths -- both work.)
EOF
