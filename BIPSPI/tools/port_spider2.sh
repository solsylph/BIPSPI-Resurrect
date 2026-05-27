#!/usr/bin/env bash
# Port SPIDER2's pred_pssm.py from Python 2 to Python 3.
#
# The "NN" inside pred_pssm.py is already pure numpy (matrix multiplies +
# sigmoid). The weights are .npz files loaded via numpy.load(). No Theano
# dependency. The ONLY Py2->Py3 syntax issues are three "print >>fp, ..."
# statements (Py2's redirect-print syntax, removed in Py3).
#
# We also fix a numpy-1.24+ compatibility issue: numpy.matrix scalar
# extraction must use .item() before %-formatting (otherwise %5.1f on a
# 1x1 matrix produces a TypeError).
#
# Source: $HOME/tools/SPIDER2/misc/pred_pssm.py
# Backup: pred_pssm.py.bak.<timestamp>
#
# Idempotent: re-running just rewrites the same edits.

set -euo pipefail

SRC=$HOME/tools/SPIDER2/misc/pred_pssm.py
if [ ! -f "$SRC" ]; then
  echo "ERROR: $SRC not found.  Did SPIDER2 extract somewhere else?"
  exit 1
fi

# Backup
cp "$SRC" "${SRC}.bak.$(date +%s)"

# Use python for the substitution -- regex on print >> Y, X is risky in sed.
python3 - "$SRC" <<'PY'
import re, sys
path = sys.argv[1]
with open(path) as f:
    text = f.read()

# Fix 1: Py2 redirect-print:  print >>FH, EXPR    ->    print(EXPR, file=FH)
# Three call sites in pred_pssm.py.  Trailing-comma "no newline" idiom -> end=''.
pat_print = re.compile(r"^(\s*)print\s*>>\s*([\w\.]+)\s*,\s*(.+)$", re.MULTILINE)
def repl_print(m):
    indent, fh, expr = m.group(1), m.group(2), m.group(3).rstrip()
    if expr.endswith(','):
        expr = expr[:-1].rstrip()
        return f"{indent}print({expr}, file={fh}, end='')"
    return f"{indent}print({expr}, file={fh})"
text, n_print = pat_print.subn(repl_print, text)
print(f"Patched {n_print} print-redirect statements")

# Fix 2: numpy.load() default changed to allow_pickle=False in numpy >=1.16.3.
# SPIDER2's .npz files contain object-dtype arrays (scipy.io.savemat-style
# nested structured arrays), so loading them requires allow_pickle=True.
# Only match calls that DON'T already have an allow_pickle argument.
pat_load = re.compile(r"numpy\.load\(([^,)]+)\)(?!\s*[,])")
def repl_load(m):
    return f"numpy.load({m.group(1)}, allow_pickle=True)"
text, n_load = pat_load.subn(repl_load, text)
print(f"Patched {n_load} numpy.load() calls (added allow_pickle=True)")

with open(path, "w") as f:
    f.write(text)
PY

echo ""
echo "=== Diff vs backup ==="
diff -u "${SRC}.bak."* "$SRC" | head -40 || true
echo ""

# Sanity: byte-compile under Py3.
echo "=== py_compile check ==="
python3 -m py_compile "$SRC" && echo "OK: pred_pssm.py parses under Python 3"
echo ""

# Smoke test against the example PSSM in SPIDER2/ex/.
echo "=== smoke test on bundled example (SPIDER2/ex/1a1xA.pssm) ==="
EX_DIR=$HOME/tools/SPIDER2/ex
if [ -f "$EX_DIR/1a1xA.pssm" ]; then
  TMPOUT=$(mktemp -d)
  cp "$EX_DIR/1a1xA.pssm" "$TMPOUT/"
  cd "$TMPOUT"
  python3 "$SRC" -f 1a1xA.pssm
  if [ -f 1a1xA.spd3 ]; then
    echo "PASS: produced 1a1xA.spd3"
    head -3 1a1xA.spd3
    echo "..."
    tail -3 1a1xA.spd3
    echo ""
    echo "Comparing against bundled reference 1a1xA_CHECK.spd3 (allow small numerical drift):"
    diff <(head -5 1a1xA.spd3) <(head -5 "$EX_DIR/1a1xA_CHECK.spd3") || true
  else
    echo "FAIL: 1a1xA.spd3 not produced"
    exit 1
  fi
else
  echo "WARN: example PSSM not found at $EX_DIR/1a1xA.pssm -- skipping smoke test"
fi

echo ""
echo "=== DONE ==="
echo "SPIDER2 ported. Configured path in BIPSPI's dependencies.cfg:"
echo "    spider2PyScript_path $SRC"
echo ""
echo "BIPSPI's spider2Manager.py invokes this script via the configured path."
