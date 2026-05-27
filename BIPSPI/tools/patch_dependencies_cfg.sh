#!/usr/bin/env bash
# Patch configFiles/cmdTool/dependencies.cfg with the cluster-specific paths.
#
# Run after Phase B completes (al2co + clustalw installed).
# Phase D (uniref90) path is also written here; if you haven't built
# the DB yet, the cfg will point at the expected location and psiblast
# will fail at runtime until Phase D is done -- that's fine, intentional.
#
# Idempotent: re-running just rewrites the values to the current ones.

set -euo pipefail

CFG=configFiles/cmdTool/dependencies.cfg
if [ ! -f "$CFG" ]; then
  echo "ERROR: $CFG not found.  Run this from ~/BIPSPI-Resurrect/BIPSPI/"
  exit 1
fi

# Resolve paths from the current environment.  Anything not found falls back
# to a bare name (works on PATH after `module load`).
PSIBLAST=$(command -v psiblast 2>/dev/null || echo "psiblast")
CDHIT=$(command -v cd-hit 2>/dev/null || echo "cd-hit")
CLUSTALW=$(command -v clustalw 2>/dev/null || echo "$HOME/miniforge3/envs/protein/bin/clustalw")
AL2CO=$HOME/tools/al2co/al2co
UNIREF90=$HOME/databases/uniref90/uniref90.fasta
SPIDER2=$HOME/tools/SPIDER2/misc/pred_pssm.py   # set after Phase C; safe placeholder

# Sanity check binaries we expect to be ready.
for bin in "$AL2CO" "$CLUSTALW"; do
  if [ ! -x "$bin" ] && ! command -v "$bin" >/dev/null 2>&1; then
    echo "WARN: $bin not found or not executable -- writing to cfg anyway"
  fi
done

echo "Patching $CFG with paths:"
echo "  psiBlastBin       $PSIBLAST"
echo "  psiBlastDB_path   $UNIREF90"
echo "  cdHitBin_path     $CDHIT"
echo "  clustalW_path     $CLUSTALW"
echo "  al2coBin_path     $AL2CO"
echo "  spider2PyScript_path $SPIDER2"

# Backup before editing.
cp "$CFG" "${CFG}.bak.$(date +%s)"

# Use python for the substitution to avoid sed escaping pain with absolute paths.
python3 - "$CFG" "$PSIBLAST" "$UNIREF90" "$CDHIT" "$CLUSTALW" "$AL2CO" "$SPIDER2" <<'PY'
import re, sys
cfg_path, psiblast, uniref90, cdhit, clustalw, al2co, spider2 = sys.argv[1:]
with open(cfg_path) as f:
    text = f.read()

# Map: cfg key -> new value.  BIPSPI's Config.py auto-strips "_path" suffix when
# reading; we write whichever key name appears in the file already.
subs = {
    "psiBlastBin":        psiblast,
    "psiBlastDB_path":    uniref90,
    "cdHitBin_path":      cdhit,
    "clustalW_path":      clustalw,
    "al2coBin_path":      al2co,
    "spider2PyScript_path": spider2,
}
for key, val in subs.items():
    # match "key  <anything>" at line start, preserve any inline comment.
    pat = re.compile(r"^(" + re.escape(key) + r")\s+\S+(.*)$", re.MULTILINE)
    new_line = r"\1 " + val + r"\2"
    text, n = pat.subn(new_line, text)
    if n == 0:
        # Key wasn't already in the file -- append it.
        text += f"\n{key} {val}\n"

with open(cfg_path, "w") as f:
    f.write(text)
print(f"Patched {cfg_path}")
PY

echo ""
echo "=== Diff vs backup ==="
diff -u "${CFG}.bak."* "$CFG" || true
echo ""
echo "Done."
