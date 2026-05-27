#!/usr/bin/env bash
# Phase B: install al2co + clustalw 1.83 from source.
#
# al2co (Pei & Grishin, swmed.edu) computes per-position conservation scores
# from a multiple-sequence alignment.  clustalw 1.83 is its alignment input.
# Both are tiny C programs with no library dependencies beyond libc.
#
# Per BIPSPI docs/repo_help.md, al2co's default name-buffer size (char[500])
# is too small for many real headers -- patch to char[1024] before building.
#
# Outputs:
#   ~/tools/al2co/al2co        (binary)
#   ~/tools/clustalw1.83/clustalw   (binary)
#
# After running this, update configFiles/cmdTool/dependencies.cfg to point at
# them (see end of script for the exact lines).

set -euo pipefail

TOOLS=$HOME/tools
mkdir -p "$TOOLS"
cd "$TOOLS"

# ---------------------------------------------------------------------------
# 1. clustalw 1.83
# ---------------------------------------------------------------------------
if [ ! -x "$TOOLS/clustalw1.83/clustalw" ]; then
  echo "=== Installing clustalw 1.83 ==="
  CW_URL="http://www.clustal.org/download/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalW/clustalw1.83.UNIX.tar.gz"
  if ! wget -q --show-progress "$CW_URL"; then
    echo "ERROR: download from $CW_URL failed. Try alternative mirrors or pull manually."
    echo "Alternative: search 'clustalw1.83.UNIX.tar.gz' and place in $TOOLS, then re-run."
    exit 1
  fi
  tar xzf clustalw1.83.UNIX.tar.gz
  cd clustalw1.83
  make
  if [ ! -x ./clustalw ]; then
    echo "ERROR: clustalw build failed -- binary not produced"
    exit 1
  fi
  echo "clustalw built at $TOOLS/clustalw1.83/clustalw"
  cd "$TOOLS"
else
  echo "clustalw already present at $TOOLS/clustalw1.83/clustalw -- skipping"
fi

# ---------------------------------------------------------------------------
# 2. al2co (with buffer patch)
# ---------------------------------------------------------------------------
if [ ! -x "$TOOLS/al2co/al2co" ]; then
  echo "=== Installing al2co (with buffer patch) ==="
  AL_URL="http://prodata.swmed.edu/download/pub/AL2CO/al2co.tar.gz"
  if ! wget -q --show-progress "$AL_URL"; then
    echo "ERROR: download from $AL_URL failed."
    echo "Alternative: pull from https://github.com search 'al2co' for mirrors."
    exit 1
  fi
  mkdir -p al2co
  tar xzf al2co.tar.gz -C al2co --strip-components=1 2>/dev/null \
    || tar xzf al2co.tar.gz -C al2co
  cd al2co

  # Patch: enlarge name buffers from 500 to 1024 (per BIPSPI docs/repo_help.md).
  # Catches the common `char xxx[500]` declarations across the .c files.
  echo "Patching name-buffer sizes 500 -> 1024..."
  shopt -s nullglob
  for f in *.c *.h; do
    [ -f "$f" ] || continue
    # be conservative: only widen *exact* char[500] declarations,
    # don't touch e.g. char[5000] or unrelated 500 constants.
    sed -i -E 's/\bchar([[:space:]]+[A-Za-z_][A-Za-z0-9_]*)\[500\]/char\1[1024]/g' "$f"
  done
  shopt -u nullglob

  make
  if [ ! -x ./al2co ]; then
    echo "ERROR: al2co build failed -- binary not produced"
    exit 1
  fi
  echo "al2co built at $TOOLS/al2co/al2co"
  cd "$TOOLS"
else
  echo "al2co already present at $TOOLS/al2co/al2co -- skipping"
fi

# ---------------------------------------------------------------------------
# 3. Print the dependencies.cfg lines you need
# ---------------------------------------------------------------------------
cat <<EOF

=== DONE ===

Now edit ~/BIPSPI-Resurrect/BIPSPI/configFiles/cmdTool/dependencies.cfg so
the following lines point at your built binaries:

    al2coBin_path     $TOOLS/al2co/al2co
    clustalW_path     $TOOLS/clustalw1.83/clustalw
    cdHitBin_path     \$(which cd-hit)
    psiBlastBin       \$(which psiblast)
    psiBlastDB_path   /path/to/uniref90.blastdb   # see Phase D download script

(cd-hit + psiblast come from \`module load\` so you can use their absolute
paths from \`which\`, or just keep them as bare \`cd-hit\` / \`psiblast\`
since they're on PATH after module load.)
EOF
