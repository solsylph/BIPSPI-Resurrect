#!/usr/bin/env python
"""Single-complex codification smoke test.

Validates that ONE complex codifies cleanly before launching a full Path-B run
(codification does NOT cache, so a failed full run wastes the whole encode).
Run from the BIPSPI root with PYTHONPATH=. so local packages import:

    PYTHONPATH=. python dbg_codify.py --wdir ~/bipspi_run/esm2_splits/wdir 4lvhbc

Expects the feature files for <prefix> to already exist under
<wdir>/computedFeatures/. Prints SUCCESS / FAILURE and exits 0 / 1.

NOTE on prefixes: BIPSPI reserves '@' as a sampling-variant tag and strips
everything after it, so prefixes here must be '@'-free (e.g. "4lvhbc", not
"4lvh@bc"). See tools/prepare_bipspi_inputs.py:make_prefix.
"""
from __future__ import absolute_import, print_function
import os
import sys
import argparse
import traceback


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="complex prefix to codify, e.g. 4lvhbc (no '@')")
    ap.add_argument("--wdir", required=True,
                    help="run wdir; computedFeatures/ lives directly under it")
    ap.add_argument("--environ", default="seq",
                    help="environType passed to OneComplexCodifier (default: seq)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if "@" in args.prefix:
        print("REFUSING: prefix %r contains '@' -- BIPSPI strips it. Use the "
              "'@'-free form (e.g. 4lvhbc)." % args.prefix, file=sys.stderr)
        return 2

    wdir = os.path.expanduser(args.wdir)
    dataRootPath = os.path.join(wdir, "computedFeatures")
    if not os.path.isdir(dataRootPath):
        print("FAILURE: no computedFeatures dir at %s" % dataRootPath, file=sys.stderr)
        return 1

    from codifyComplexes.codifyOneComplex import OneComplexCodifier

    cod = OneComplexCodifier(dataRootPath=dataRootPath, environType=args.environ,
                             verbose=args.verbose)
    try:
        result = cod.codifyComplex(args.prefix)
    except Exception:
        print("FAILURE: codifyComplex(%r) raised:" % args.prefix, file=sys.stderr)
        traceback.print_exc()
        return 1

    if result is None:
        print("FAILURE: codifyComplex(%r) returned None" % args.prefix, file=sys.stderr)
        return 1

    # ComplexCodified exposes the encoded pair table; report its shape if present.
    shape = getattr(getattr(result, "pairsDf", None), "shape", None)
    print("SUCCESS: codified %s%s" % (args.prefix,
          (" (pairsDf shape %s)" % (shape,)) if shape is not None else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
