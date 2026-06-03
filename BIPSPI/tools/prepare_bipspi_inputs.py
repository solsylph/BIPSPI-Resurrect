#!/usr/bin/env python3
"""Adapt ESMFold-multimer pipeline outputs into BIPSPI training inputs.

Reads:
  - splits dir with train.json / val.json / test.json (from
    code/esmfold_multimer/scripts/data/05_create_splits.py)
  - structures dir with {PDB_ID}_assembly1.cif (from 03_download_structures.py)

Writes:
  - <output-dir>/pdbs/{PREFIX}_l_b.pdb, {PREFIX}_l_u.pdb,
                       {PREFIX}_r_b.pdb, {PREFIX}_r_u.pdb
        one set per complex; _u files are symlinks to _b (no unbound
        coords in our pipeline), as documented in docs/repo_help.md.
  - <output-dir>/folds.json
        single-fold JSON consumable via BIPSPI's `--N_KFOLD <foldsFile>` flag:
        [{"train": [<train+val prefixes>], "test": [<test prefixes>]}]
  - <output-dir>/manifest.json
        map prefix -> (pdb_id, chain_a, chain_b, split) for traceability.

Prefix convention (per BIPSPI docs/repo_help.md):
  - UPPERCASE prefix  -> complex is evaluated (SKIP_LOWER_PREDICTION
    in trainAndTest.py drops lowercase prefixes at the last step).
  - lowercase prefix  -> train-only complex.
By default: train + val complexes are lowercased; test complexes are uppercased.
That gives a single-run baseline producing test-set metrics directly comparable
to 07b's RF baseline.

Disambiguation: (pdb_id, chain_a, chain_b) -> "{PDB}@{chainA}{chainB}".
BIPSPI uses prefix[:4].isupper() for the eval gate so the @-suffix doesn't
break case detection.

Multi-copy chains: if the biological assembly CIF contains multiple physical
copies of chain_a (e.g. C2 homo-dimer where label "A" appears twice), all
copies are kept in the ligand PDB. BIPSPI's struct/seq feature pipeline treats
this as a single partner with repeated sequence.

Dependencies: gemmi, tqdm (both standard in the ESMFold pipeline env).
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import gemmi
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def make_prefix(pdb_id: str, chain_a: str, chain_b: str, evaluated: bool) -> str:
    """Build a BIPSPI-style prefix. Uppercase if evaluated, lowercase otherwise.

    The 4-char PDB ID slice is the part BIPSPI's evaluation gate checks
    (prefix[:4].isupper()).  Chain ids are appended after '@' for uniqueness.
    """
    pdb = pdb_id.upper() if evaluated else pdb_id.lower()
    chains = (chain_a + chain_b)
    chains = chains.upper() if evaluated else chains.lower()
    return f"{pdb}@{chains}"


def extract_chain_to_pdb(
    cif_path: Path,
    chain_id: str,
    out_pdb_path: Path,
) -> bool:
    """Extract a single chain (matching auth_asym_id) from an mmCIF assembly
    file and write it as a PDB.  Returns True on success, False if no
    matching chain was found.
    """
    structure = gemmi.read_structure(str(cif_path))
    if len(structure) == 0:
        return False
    # gemmi reads each model; we use model 0 (assembly1 has only one model).
    src_model = structure[0]

    # Build new structure containing only chains whose name (auth_asym_id under
    # the gemmi default read) matches chain_id.  Multiple physical copies all
    # share the same name and are all kept.
    new_struct = gemmi.Structure()
    new_struct.cell = structure.cell
    new_struct.spacegroup_hm = structure.spacegroup_hm
    # gemmi >=0.6 renamed Model.name -> Model.num (int); Model() constructor
    # still takes a string identifier.
    new_model = gemmi.Model(str(src_model.num))

    n_copies = 0
    for chain in src_model:
        if chain.name == chain_id:
            new_model.add_chain(chain.clone())
            n_copies += 1

    if n_copies == 0:
        return False

    new_struct.add_model(new_model)
    out_pdb_path.parent.mkdir(parents=True, exist_ok=True)
    new_struct.write_pdb(str(out_pdb_path))
    return True


def symlink_or_copy(src: Path, dst: Path) -> None:
    """Make dst point at src.  Prefer hard-link/symlink; fall back to copy
    (Windows without admin doesn't allow symlinks, Linux/cluster is fine).
    """
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.symlink(src.name, dst)  # relative symlink (same dir)
        return
    except (OSError, NotImplementedError):
        pass
    try:
        os.link(src, dst)
        return
    except (OSError, NotImplementedError):
        pass
    import shutil
    shutil.copyfile(src, dst)


def process_complex(
    cx: Dict,
    structures_dir: Path,
    out_pdbs_dir: Path,
    evaluated: bool,
) -> Optional[str]:
    """Process one complex.  Returns prefix on success, None on skip."""
    pdb_id = cx["pdb_id"]
    chain_a = cx["chain_a_id"]
    chain_b = cx["chain_b_id"]
    prefix = make_prefix(pdb_id, chain_a, chain_b, evaluated)

    cif_path = structures_dir / f"{pdb_id}_assembly1.cif"
    if not cif_path.exists():
        cif_path = structures_dir / f"{pdb_id.lower()}_assembly1.cif"
    if not cif_path.exists():
        logger.warning(f"{prefix}: missing CIF {pdb_id}_assembly1.cif")
        return None

    lig_b = out_pdbs_dir / f"{prefix}_l_b.pdb"
    rec_b = out_pdbs_dir / f"{prefix}_r_b.pdb"
    lig_u = out_pdbs_dir / f"{prefix}_l_u.pdb"
    rec_u = out_pdbs_dir / f"{prefix}_r_u.pdb"

    if all(p.exists() for p in (lig_b, rec_b, lig_u, rec_u)):
        return prefix  # already prepared, idempotent

    ok_a = extract_chain_to_pdb(cif_path, chain_a, lig_b)
    if not ok_a:
        logger.warning(f"{prefix}: chain_a '{chain_a}' not found in {cif_path.name}")
        return None
    ok_b = extract_chain_to_pdb(cif_path, chain_b, rec_b)
    if not ok_b:
        logger.warning(f"{prefix}: chain_b '{chain_b}' not found in {cif_path.name}")
        # roll back ligand
        if lig_b.exists():
            lig_b.unlink()
        return None

    symlink_or_copy(lig_b, lig_u)
    symlink_or_copy(rec_b, rec_u)
    return prefix


def load_split(path: Path) -> List[Dict]:
    with open(path) as fh:
        data = json.load(fh)
    return data["complexes"]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--splits-dir", required=True, type=Path,
                        help="Path to directory with train.json/val.json/test.json")
    parser.add_argument("--structures-dir", required=True, type=Path,
                        help="Path to directory with *_assembly1.cif files")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Where to write BIPSPI inputs (will create pdbs/, folds.json, manifest.json)")
    parser.add_argument("--evaluate", choices=["test", "val+test", "all"], default="test",
                        help="Which splits get UPPERCASE prefixes (evaluated). "
                             "Default 'test' = train+val are train-only, test is evaluated.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only first N complexes per split (smoke testing).")
    args = parser.parse_args()

    splits_dir = args.splits_dir.expanduser().resolve()
    structures_dir = args.structures_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    pdbs_dir = output_dir / "pdbs"
    pdbs_dir.mkdir(parents=True, exist_ok=True)

    splits = {}
    for name in ("train", "val", "test"):
        p = splits_dir / f"{name}.json"
        if not p.exists():
            sys.exit(f"Missing split file: {p}")
        splits[name] = load_split(p)
        if args.limit:
            splits[name] = splits[name][: args.limit]
        logger.info(f"Loaded {name}: {len(splits[name])} complexes")

    eval_splits = {
        "test":     {"test"},
        "val+test": {"val", "test"},
        "all":      {"train", "val", "test"},
    }[args.evaluate]

    manifest: Dict[str, Dict] = {}
    fold_train: List[str] = []
    fold_test: List[str] = []

    for split_name, complexes in splits.items():
        evaluated = split_name in eval_splits
        for cx in tqdm(complexes, desc=f"prepare {split_name}", unit="cx"):
            prefix = process_complex(cx, structures_dir, pdbs_dir, evaluated)
            if prefix is None:
                continue
            manifest[prefix] = {
                "pdb_id":  cx["pdb_id"],
                "chain_a": cx["chain_a_id"],
                "chain_b": cx["chain_b_id"],
                "split":   split_name,
                "evaluated": evaluated,
            }
            if split_name == "test":
                fold_test.append(prefix)
                fold_train.append(prefix)  # test must appear in train for BIPSPI's CV plumbing; isLastStep skips lowercase
            else:
                fold_train.append(prefix)

    # Write folds JSON consumable by `--N_KFOLD <foldsFile>` in BIPSPI.
    # Single fold: train = train+val (lowercase) + test (uppercase, included for
    # BIPSPI's cross-validation plumbing which expects test prefixes to also
    # appear in the train list); test = uppercase test prefixes.
    # BIPSPI's SKIP_LOWER_PREDICTION + last-step gate produces metrics on
    # uppercase prefixes only -> our test set.
    folds = [{"train": fold_train, "test": fold_test}]
    with open(output_dir / "folds.json", "w") as fh:
        json.dump(folds, fh, indent=2)
    with open(output_dir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)

    logger.info(f"Wrote {len(manifest)} complexes to {pdbs_dir}")
    logger.info(f"Wrote folds.json: train={len(fold_train)} prefixes, test={len(fold_test)} prefixes")
    logger.info(f"Wrote manifest.json: {output_dir / 'manifest.json'}")


if __name__ == "__main__":
    main()
