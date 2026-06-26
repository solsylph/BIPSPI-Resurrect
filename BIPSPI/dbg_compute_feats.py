#!/usr/bin/env python
"""Run ONLY the feature-compute step, then stop (no codify / train / eval).

Used to regenerate feature files in place -- e.g. after fixing the ESM2 resId
bug, delete the stale esm2/*.esm2.tab.gz and re-run this to rewrite them with
correct PDB residue ids.  All already-cached features (contact maps, seq->struct
maps, ...) are skipped via their os.path.isfile checks, so only the deleted
files are recomputed.  This lets us validate one complex with dbg_codify BEFORE
committing to the full (non-cached) codify pass.

Mirrors the argument setup of generateBIPSPIModel.py's __main__ exactly, but
calls computeFeatures() instead of main().

    PYTHONPATH=. python dbg_compute_feats.py \
        --modelType seq \
        --pdbsIndir ~/bipspi_run/esm2_splits/pdbs \
        --N_KFOLD ~/bipspi_run/esm2_splits/folds.json \
        --wdir ~/bipspi_run/esm2_splits/wdir \
        --ncpu 16
"""
from __future__ import absolute_import, print_function

import generateBIPSPIModel as g
from Config import Configuration

if __name__ == "__main__":
    parser = Configuration.getArgParser()
    parser.modify_field("pdbsIndir", help="Directory where training pdbs are located", _type=Configuration.file_path)
    parser.modify_field("wdir", help="Directory where partial results and final results will be saved", _type=Configuration.file_path)
    parser.modify_field("tmp", help="Temporary directory", _type=Configuration.file_path)
    parser.modify_field("ncpu", help="Number of cpus for feature computation")
    parser.modify_field("modelType", help="struct, seq or mixed", choices=["struct", "mixed", "seq"])
    parser.modify_field("N_KFOLD", help="Cross-validation type. -1 leave-one-out, positive k, or path to folds.json", _type=Configuration.int_or_filePath)
    parser.modify_field("scopeFamiliesFname", help="Filename containing the families of the protein chains", _type=Configuration.file_path)

    parser.parse_args()

    g.computeFeatures()
    print("FEATURE COMPUTE DONE (codify/train/eval intentionally skipped)")
