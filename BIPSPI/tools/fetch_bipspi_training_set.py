#!/usr/bin/env python3
"""Fetch + prepare the canonical BIPSPI training set from info_HEDt.tab.

Two data sources:
  1. RCSB Protein Data Bank — for the ~2403 lowercase (train-only augmentation)
     entries.  Fetches each PDB, extracts the chains specified in info_HEDt.tab,
     produces BIPSPI-format PDBID_l_u.pdb / PDBID_r_u.pdb.  Symlinks
     _b -> _u per BIPSPI docs convention (no separate bound conformation
     available for these entries).

  2. Protein Docking Benchmark 5 (https://zlab.wenglab.org/benchmark/) — for the
     ~228 uppercase (evaluated) entries.  These come pre-prepared as 4-file
     BIPSPI-format bundles (PDBID_l_b.pdb, PDBID_l_u.pdb, PDBID_r_b.pdb,
     PDBID_r_u.pdb -- BIPSPI's naming convention was modelled on Benchmark 5).
     You download Benchmark 5 manually and point --benchmark5-dir at the
     extracted directory.

Outputs (all under --out-dir):
  pdbs/                      -- the BIPSPI pdbsIndir (PDBID_{l,r}_{b,u}.pdb per complex)
  info_HEDt.noheader.tab     -- header-stripped scopes file for --scopeFamiliesFname
  manifest.json              -- per-complex status (source, ok, message)
  failed.txt                 -- list of complexes that couldn't be fetched/extracted
  raw/                       -- cached raw RCSB PDB downloads (can be deleted after)

Usage (smoke test with 5 of each kind):
  python tools/fetch_bipspi_training_set.py \
      --info-tab docs/info_HEDt.tab \
      --benchmark5-dir ~/benchmark5 \
      --out-dir ~/bipspi_run/canonical_smoke \
      --limit 5

Usage (full set):
  python tools/fetch_bipspi_training_set.py \
      --info-tab docs/info_HEDt.tab \
      --benchmark5-dir ~/benchmark5 \
      --out-dir ~/bipspi_run/canonical \
      --workers 8

Resumable + idempotent: re-running skips complexes already present in pdbs/.
"""

import argparse
import json
import logging
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import gemmi
import requests
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


RCSB_URL_PDB = "https://files.rcsb.org/download/{pdbid_upper}.pdb"


def parse_info_tab(path: Path) -> List[Dict]:
    """Parse info_HEDt.tab into a list of complex records.
    Tolerates a header row ('pdbId\\tchainIds1\\t...') if present.
    Returns: list of {pdb_id, chain_l, chain_r, scopes_l, scopes_r, evaluated}.
    """
    entries = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            pdb_id = parts[0].strip()
            if pdb_id == "pdbId":  # header row
                continue
            entries.append({
                "pdb_id":    pdb_id,
                "chain_l":   parts[1].strip(),
                "chain_r":   parts[2].strip(),
                "scopes_l":  parts[3].strip() if len(parts) > 3 else "",
                "scopes_r":  parts[4].strip() if len(parts) > 4 else "",
                "evaluated": pdb_id.isupper(),
            })
    return entries


def fetch_pdb_from_rcsb(pdb_id: str, dest: Path, session: requests.Session,
                        timeout: int = 30) -> bool:
    """Download a raw PDB file from RCSB.  Idempotent.  Returns True on success."""
    if dest.exists() and dest.stat().st_size > 0:
        return True
    url = RCSB_URL_PDB.format(pdbid_upper=pdb_id.upper())
    try:
        r = session.get(url, timeout=timeout)
        if r.status_code != 200:
            return False
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(r.content)
        return True
    except requests.RequestException:
        return False


def extract_chains_to_pdb(src_pdb: Path, chain_spec: str, dest_pdb: Path) -> bool:
    """Extract one or more chains from src_pdb into dest_pdb.

    chain_spec is a concatenated chain id string (e.g. 'A', 'AB').
    '+' means the no-id (blank) chain per info_HEDt.tab convention.

    Returns True if at least one matching chain was found and written.
    """
    if dest_pdb.exists() and dest_pdb.stat().st_size > 0:
        return True

    # Each char is one chain id; '+' -> match blank chain name.
    target_ids = set()
    match_blank = False
    for c in chain_spec:
        if c == "+":
            match_blank = True
        else:
            target_ids.add(c)

    try:
        structure = gemmi.read_structure(str(src_pdb))
    except RuntimeError:
        return False
    if len(structure) == 0:
        return False
    src_model = structure[0]

    new_struct = gemmi.Structure()
    new_struct.cell = structure.cell
    new_struct.spacegroup_hm = structure.spacegroup_hm
    # gemmi >=0.6 renamed Model.name -> Model.num (int); Model() constructor
    # still takes a string identifier.
    new_model = gemmi.Model(str(src_model.num))

    n_chains = 0
    for chain in src_model:
        is_blank = chain.name.strip() == ""
        if chain.name in target_ids or (match_blank and is_blank):
            new_model.add_chain(chain.clone())
            n_chains += 1

    if n_chains == 0:
        return False

    new_struct.add_model(new_model)
    dest_pdb.parent.mkdir(parents=True, exist_ok=True)
    new_struct.write_pdb(str(dest_pdb))
    return True


def link_or_copy(src: Path, dst: Path) -> None:
    """Make dst point at src.  Uses relative symlink when in the same dir,
    absolute symlink when cross-dir, falls back to copy if symlinks aren't
    supported (e.g. Windows without dev mode)."""
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        if src.parent.resolve() == dst.parent.resolve():
            os.symlink(src.name, dst)  # relative
        else:
            os.symlink(str(src.resolve()), dst)  # absolute
        return
    except (OSError, NotImplementedError):
        pass
    shutil.copyfile(src, dst)


def process_lowercase_entry(entry: Dict, raw_pdb_dir: Path, out_pdbs_dir: Path,
                            session: requests.Session) -> Tuple[str, bool, str]:
    """RCSB fetch + chain extraction + _b -> _u symlinks for one train-only complex.
    Returns (prefix, success, message)."""
    prefix = entry["pdb_id"]  # case preserved (lowercase)
    raw_pdb = raw_pdb_dir / f"{prefix}.pdb"

    lig_u = out_pdbs_dir / f"{prefix}_l_u.pdb"
    rec_u = out_pdbs_dir / f"{prefix}_r_u.pdb"
    lig_b = out_pdbs_dir / f"{prefix}_l_b.pdb"
    rec_b = out_pdbs_dir / f"{prefix}_r_b.pdb"

    if all(p.exists() for p in (lig_u, rec_u, lig_b, rec_b)):
        return (prefix, True, "already prepared")

    if not fetch_pdb_from_rcsb(prefix, raw_pdb, session):
        return (prefix, False, "rcsb fetch failed")

    if not extract_chains_to_pdb(raw_pdb, entry["chain_l"], lig_u):
        return (prefix, False, f"chain_l '{entry['chain_l']}' not found in {prefix}.pdb")
    if not extract_chains_to_pdb(raw_pdb, entry["chain_r"], rec_u):
        if lig_u.exists():
            lig_u.unlink()
        return (prefix, False, f"chain_r '{entry['chain_r']}' not found in {prefix}.pdb")

    link_or_copy(lig_u, lig_b)
    link_or_copy(rec_u, rec_b)
    return (prefix, True, "ok")


def process_uppercase_entry(entry: Dict, benchmark5_dir: Path,
                            out_pdbs_dir: Path) -> Tuple[str, bool, str]:
    """Locate the 4-file BIPSPI-format bundle for one Benchmark 5 complex and
    link/copy into the output pdbsIndir.  Returns (prefix, success, message)."""
    prefix = entry["pdb_id"]  # uppercase preserved

    # Benchmark 5 uses BIPSPI's exact naming.
    sources = {suffix: benchmark5_dir / f"{prefix}_{suffix}.pdb"
               for suffix in ("l_b", "l_u", "r_b", "r_u")}
    missing = [k for k, p in sources.items() if not p.exists()]
    if missing:
        return (prefix, False, f"benchmark5 files missing: {missing}")

    for suffix, src in sources.items():
        dst = out_pdbs_dir / f"{prefix}_{suffix}.pdb"
        if dst.exists():
            continue
        link_or_copy(src, dst)
    return (prefix, True, "ok")


def write_noheader_scopes(entries: List[Dict], out_path: Path) -> None:
    """Write info_HEDt.tab without the header, suitable for BIPSPI's
    --scopeFamiliesFname.  Tab-separated, 5 columns."""
    with open(out_path, "w") as fh:
        for e in entries:
            fh.write(f"{e['pdb_id']}\t{e['chain_l']}\t{e['chain_r']}\t"
                     f"{e['scopes_l']}\t{e['scopes_r']}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--info-tab", required=True, type=Path,
                        help="Path to info_HEDt.tab (with or without header)")
    parser.add_argument("--benchmark5-dir", type=Path, default=None,
                        help="Directory containing extracted Benchmark 5 files "
                             "(PDBID_l_b.pdb etc.). Required to process uppercase entries.")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Where to write BIPSPI inputs (pdbs/, manifest.json, etc.)")
    parser.add_argument("--raw-pdb-dir", type=Path, default=None,
                        help="Where to cache RCSB downloads (default: <out-dir>/raw)")
    parser.add_argument("--workers", type=int, default=8,
                        help="Parallel HTTP workers for RCSB fetch (default 8)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only first N complexes per source (smoke testing)")
    parser.add_argument("--skip-uppercase", action="store_true",
                        help="Don't process uppercase / Benchmark 5 entries")
    parser.add_argument("--skip-lowercase", action="store_true",
                        help="Don't process lowercase / RCSB entries")
    args = parser.parse_args()

    info_tab = args.info_tab.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    raw_pdb_dir = (args.raw_pdb_dir or (out_dir / "raw")).expanduser().resolve()
    pdbs_dir = out_dir / "pdbs"
    pdbs_dir.mkdir(parents=True, exist_ok=True)
    raw_pdb_dir.mkdir(parents=True, exist_ok=True)

    entries = parse_info_tab(info_tab)
    logger.info(f"Parsed {len(entries)} entries from {info_tab.name}")

    scopes_out = out_dir / "info_HEDt.noheader.tab"
    write_noheader_scopes(entries, scopes_out)
    logger.info(f"Wrote header-stripped scopes file: {scopes_out}")

    lowercase = [e for e in entries if not e["evaluated"]]
    uppercase = [e for e in entries if e["evaluated"]]
    if args.limit:
        lowercase = lowercase[: args.limit]
        uppercase = uppercase[: args.limit]
    logger.info(f"To process: {len(lowercase)} lowercase (RCSB), "
                f"{len(uppercase)} uppercase (Benchmark 5)")

    manifest = {}
    failed = []

    # Uppercase (Benchmark 5) first -- local I/O only, fast.
    if not args.skip_uppercase:
        if not args.benchmark5_dir:
            logger.warning("--benchmark5-dir not provided; skipping uppercase entries")
        else:
            b5dir = args.benchmark5_dir.expanduser().resolve()
            if not b5dir.exists():
                logger.error(f"--benchmark5-dir does not exist: {b5dir}")
                sys.exit(2)
            for e in tqdm(uppercase, desc="benchmark5", unit="cx"):
                prefix, ok, msg = process_uppercase_entry(e, b5dir, pdbs_dir)
                manifest[prefix] = {
                    "source":    "benchmark5",
                    "ok":        ok,
                    "message":   msg,
                    "chain_l":   e["chain_l"],
                    "chain_r":   e["chain_r"],
                    "evaluated": True,
                }
                if not ok:
                    failed.append((prefix, msg))

    # Lowercase (RCSB) -- network-bound, parallelize.
    if not args.skip_lowercase and lowercase:
        with requests.Session() as session:
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                futures = {
                    pool.submit(process_lowercase_entry, e, raw_pdb_dir, pdbs_dir,
                                session): e
                    for e in lowercase
                }
                for fut in tqdm(as_completed(futures), total=len(futures),
                                desc="rcsb", unit="cx"):
                    e = futures[fut]
                    try:
                        prefix, ok, msg = fut.result()
                    except Exception as exc:
                        prefix, ok, msg = e["pdb_id"], False, f"exception: {exc}"
                    manifest[prefix] = {
                        "source":    "rcsb",
                        "ok":        ok,
                        "message":   msg,
                        "chain_l":   e["chain_l"],
                        "chain_r":   e["chain_r"],
                        "evaluated": False,
                    }
                    if not ok:
                        failed.append((prefix, msg))

    with open(out_dir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    if failed:
        with open(out_dir / "failed.txt", "w") as fh:
            for prefix, msg in sorted(failed):
                fh.write(f"{prefix}\t{msg}\n")

    n_ok = sum(1 for v in manifest.values() if v["ok"])
    logger.info(f"Done: {n_ok}/{len(manifest)} prepared, {len(failed)} failed")
    if failed:
        logger.info(f"Failures listed in {out_dir / 'failed.txt'}")
    logger.info(f"BIPSPI pdbsIndir ready at: {pdbs_dir}")
    logger.info(f"BIPSPI scopesFamilies file ready at: {scopes_out}")


if __name__ == "__main__":
    main()
