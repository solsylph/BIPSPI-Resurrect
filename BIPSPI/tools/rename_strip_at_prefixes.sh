#!/usr/bin/env bash
# Strip the reserved '@' character out of an existing Path-B run's artifacts so
# they match the '@'-free prefixes that the patched prepare_bipspi_inputs.py now
# emits (e.g. "4lvh@bc" -> "4lvhbc").  This lets us keep the ~1hr of computed
# ESM2 features (rename in place) instead of recomputing.
#
# BIPSPI reserves '@' as a "same-complex sampling variant" tag and strips
# everything after it in ~15 places, which breaks feature-file lookup and
# collapses distinct chain-pairs of one pdb.  See STATUS.md s15.
#
# The transform is purely "delete the '@' character": '@' never appears in these
# names/contents except as the prefix separator, so stripping it is exactly the
# new scheme.  Bijective and collision-free (4-char PDB-id boundary).
#
# Usage:
#   bash tools/rename_strip_at_prefixes.sh <run_dir>            # DRY RUN (default)
#   bash tools/rename_strip_at_prefixes.sh <run_dir> --apply    # actually do it
#
# <run_dir> is e.g. ~/bipspi_run/esm2_splits  (contains pdbs/, wdir/, folds.json)
set -euo pipefail

RUN_DIR="${1:-}"
APPLY=0
[[ "${2:-}" == "--apply" ]] && APPLY=1

if [[ -z "$RUN_DIR" || ! -d "$RUN_DIR" ]]; then
  echo "usage: $0 <run_dir> [--apply]   (run_dir must exist)" >&2
  exit 2
fi
RUN_DIR="$(cd "$RUN_DIR" && pwd)"

WDIR="$RUN_DIR/wdir"
say() { if [[ $APPLY -eq 1 ]]; then echo "[apply] $*"; else echo "[dry  ] $*"; fi; }

echo "== run_dir: $RUN_DIR  (apply=$APPLY) =="

# ---------------------------------------------------------------------------
# 1. Rename files and directories whose basename contains '@'.
#    -depth processes children before parents so renaming a dir doesn't
#    invalidate paths to files still inside it.
# ---------------------------------------------------------------------------
renamed=0
skipped=0
for base in "$RUN_DIR/pdbs" "$WDIR/computedFeatures"; do
  [[ -d "$base" ]] || { echo "  (no $base, skipping)"; continue; }
  while IFS= read -r -d '' path; do
    dir="$(dirname "$path")"
    name="$(basename "$path")"
    newname="${name//@/}"
    newpath="$dir/$newname"
    if [[ "$newpath" == "$path" ]]; then
      continue
    fi
    if [[ -e "$newpath" ]]; then
      echo "  !! TARGET EXISTS, refusing to clobber: $newpath  (from $path)" >&2
      skipped=$((skipped+1))
      continue
    fi
    say "mv '$path' -> '$newpath'"
    if [[ $APPLY -eq 1 ]]; then mv "$path" "$newpath"; fi
    renamed=$((renamed+1))
  done < <(find "$base" -depth -name '*@*' -print0)
done
echo "  files/dirs to rename: $renamed   refused(clobber): $skipped"

# ---------------------------------------------------------------------------
# 1b. Re-point symlinks whose stored TARGET still contains '@'.  `mv` renames a
#     symlink file but NOT its target string.  The _u.pdb files are relative
#     symlinks to the _b.pdb files (prepare_bipspi_inputs.py: os.symlink(src.name,
#     dst)), so after renaming both, 4lvhbc_l_u.pdb still points at the OLD
#     4lvh@bc_l_b.pdb -> dangling -> FileNotFoundError at feature-compute time.
#     Strip '@' from the target so it points at the renamed _b file.
# ---------------------------------------------------------------------------
relinked=0
for base in "$RUN_DIR/pdbs" "$WDIR/computedFeatures"; do
  [[ -d "$base" ]] || continue
  while IFS= read -r -d '' link; do
    tgt="$(readlink "$link")"
    case "$tgt" in
      *@*)
        newtgt="${tgt//@/}"
        say "relink '$link' -> '$newtgt' (was '$tgt')"
        if [[ $APPLY -eq 1 ]]; then ln -sfn "$newtgt" "$link"; fi
        relinked=$((relinked+1))
        ;;
    esac
  done < <(find "$base" -type l -print0)
done
echo "  symlinks re-pointed: $relinked"

# ---------------------------------------------------------------------------
# 2. Rewrite JSON prefix lists (folds.json, manifest.json).  '@' only ever
#    appears inside prefix strings, so a literal global delete is safe.
# ---------------------------------------------------------------------------
for jf in "$RUN_DIR/folds.json" "$RUN_DIR/manifest.json"; do
  [[ -f "$jf" ]] || continue
  if grep -q '@' "$jf"; then
    say "strip '@' in $(basename "$jf") (backup -> $(basename "$jf").bak)"
    if [[ $APPLY -eq 1 ]]; then
      cp "$jf" "$jf.bak"
      sed -i 's/@//g' "$jf"
    fi
  else
    echo "  $(basename "$jf"): no '@', nothing to do"
  fi
done

# ---------------------------------------------------------------------------
# 3. Delete stale codified pickles (regenerated on next run; codify is not
#    cached anyway).
# ---------------------------------------------------------------------------
if [[ -d "$WDIR" ]]; then
  while IFS= read -r -d '' pkl; do
    say "rm stale $pkl"
    if [[ $APPLY -eq 1 ]]; then rm -f "$pkl"; fi
  done < <(find "$WDIR" \( -name '*.train.pkl.gz' -o -name '*.predict.pkl.gz' \) -print0)
fi

echo
if [[ $APPLY -eq 1 ]]; then
  echo "== DONE (applied). Next: validate one complex =="
  echo "   PYTHONPATH=. python dbg_codify.py --wdir $WDIR <prefix>   # expect SUCCESS"
else
  echo "== DRY RUN only. Re-run with --apply to make changes. =="
fi
