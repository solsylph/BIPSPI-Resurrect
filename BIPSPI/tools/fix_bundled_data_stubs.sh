#!/usr/bin/env bash
# Replace 12-byte placeholder *_u.pdb stubs with proper symlinks to *_b.pdb.
#
# BIPSPI's docs/trainingPDBsExample (and some real Protein Docking Benchmark
# entries) ship "unbound" PDB files as 12-byte text stubs containing just the
# target bound filename, instead of real symlinks.  Biopython's PDBParser sees
# zero atoms and returns an empty Structure, which crashes downstream contact
# map / bound-unbound mapping code.
#
# Per BIPSPI's documented convention (docs/repo_help.md):
#   "if no unbound pdbs available, use symlinks
#    prefix_r_b.pdb --> prefix_r_u.pdb
#    prefix_l_b.pdb --> prefix_l_u.pdb"
# This script implements that convention idempotently for any directory.
#
# Usage:
#   bash tools/fix_bundled_data_stubs.sh [DIR]
# Default DIR is ./docs/trainingPDBsExample
#
# Idempotent: re-running only touches files that are still stubs.

set -euo pipefail

DIR="${1:-./docs/trainingPDBsExample}"

if [ ! -d "$DIR" ]; then
  echo "ERROR: directory not found: $DIR" >&2
  exit 1
fi

echo "Scanning $DIR for stub _u.pdb files (<100 bytes, not already symlinks)..."
fixed=0
skipped=0
for stub in $(find "$DIR" -maxdepth 1 -name "*_u.pdb" -size -100c -type f 2>/dev/null); do
  target_b="${stub/_u.pdb/_b.pdb}"
  if [ ! -s "$target_b" ]; then
    echo "  SKIP: $stub -- no bound counterpart $target_b" >&2
    skipped=$((skipped+1))
    continue
  fi
  echo "  fixing $stub -> $(basename "$target_b")"
  rm "$stub"
  ln -s "$(basename "$target_b")" "$stub"
  fixed=$((fixed+1))
done

echo ""
echo "Done. Fixed: $fixed.  Skipped (no _b counterpart): $skipped."

# Suggest cache cleanup for the complexes we just touched.
if [ "$fixed" -gt 0 ]; then
  echo ""
  echo "If you have a previous wdir at /tmp/test_bipspi_seq_v2 (or any --wdir"
  echo "you reuse), you may want to clear cached contact maps for the"
  echo "fixed complexes so they get recomputed against the now-real coords:"
  echo ""
  echo "  WDIR=/tmp/test_bipspi_seq_v2"
  echo "  for f in \$(find $DIR -name '*_u.pdb' -type l); do"
  echo "      prefix=\$(basename \"\$f\" | cut -d_ -f1)"
  echo "      rm -rf \$WDIR/computedFeatures/common/contactMaps/\${prefix}_*"
  echo "  done"
fi
