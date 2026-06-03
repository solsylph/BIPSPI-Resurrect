# BIPSPI Python 3 Port — Status & Cluster Handoff

**Upstream:** [rsanchezgarc/BIPSPI](https://github.com/rsanchezgarc/BIPSPI) (Python 2.7)
**Our fork:** [solsylph/BIPSPI-Resurrect](https://github.com/solsylph/BIPSPI-Resurrect)
**Purpose:** sequence-mode resurrection on the ESMFold-multimer mmseqs2-clustered splits, as a non-RF baseline after the RF-on-ESM2 architectural ceiling was confirmed (pair-AUROC ~0.55 on test).
**Last update:** 2026-06-03 (after Phase D interrupt + restart, canonical fetch complete, three rounds of runtime fixes landed).

> **Repo layout:** GitHub repo has `BIPSPI/` as a subdirectory at its root (the local git repo was reinitialised at `E:\BIPSPI-Resurrect\` instead of `E:\BIPSPI-Resurrect\BIPSPI\`). On the cluster the code lives at `~/BIPSPI-Resurrect/BIPSPI/`. All commands assume `cd ~/BIPSPI-Resurrect/BIPSPI` first.

---

## Progress at a glance

| Phase | Subject | Status |
|---|---|---|
| Port | Py2.7 → Py3 mechanical port (21 files) | ✅ |
| Critic | Two P0 regressions found + fixed (alignment truncation, empty-input crash) | ✅ |
| A | Extend `protein` conda env + module-load BLAST+/CD-HIT | ✅ |
| B | Build al2co + install clustalw 2.1 (via bioconda) | ✅ |
| Cfg patch | Update `dependencies.cfg` with cluster paths | ✅ |
| C | Port SPIDER2 to Py3.10 (numpy-only, option a-i) | ✅ |
| D | Download uniref90 + makeblastdb | ⏳ **restarting** (interrupted 2026-05-27; restarted 2026-06-03 ~13:11) |
| Runtime fixes | Biopython 1.80 removals + Py2/Py3 gzip mode + gemmi 0.6 API | ✅ (3 rounds applied) |
| Canonical data | Fetch BIPSPI's published training set (info_HEDt.tab + Benchmark 5 + RCSB) | ✅ (2610/2611 prepared) |
| Smoke seq | `--modelType seq` against `./docs/trainingPDBsExample` | ⏳ blocked on Phase D |
| Canonical run | Full 2611-complex BIPSPI seq training | blocked on smoke |
| Re-baseline | Our splits + 07b-comparable evaluation | blocked on canonical run |

---

## 1. Mechanical Python 2 → 3 port (complete)

21 files edited under `git diff HEAD` (+83/-63 lines) against the upstream Py2.7 clone. Translations:

| Hazard | Files | Replacement |
|---|---|---|
| `Bio.pairwise2` + `Bio.SubsMat.MatrixInfo` | 4 | `Bio.Align.PairwiseAligner(mode="local")` + `substitution_matrices.load("BLOSUM62")` |
| `dict.iteritems()` | 4 | `dict.items()` |
| `reduce(...)` | 1 | `from functools import reduce` |
| `isinstance(x, file)` | 1 (2 hits) | `isinstance(x, io.IOBase)` |
| `.next()` method | 1 | `next(iter)` |
| `import StringIO` | 1 | `from io import StringIO` |
| `pd.set_option('precision', N)` | 4 | `pd.set_option('display.precision', N)` |
| `df.ix[]` + `df.append(..., ignore_index=True)` | 8 | `df.iloc[]` + `pd.concat([df, df.iloc[[-1]]], ignore_index=True)` |

`python -m compileall .` passes for 122 of 123 files.  The 1 failure (`monitorScreenlog.py` line 28) is a **pre-existing source bug** that broke in Python 2.7 too — out of scope.

## 2. Critic/implementor pair-programming round (complete)

After the mechanical port, a critic review found two correctness regressions:

- **P0-1 (silent wrong residue mappings):** `PairwiseAligner` in local mode returns the aligned strings TRUNCATED to the local match region; `pairwise2.align.localds` returned the full input sequences with `-` padding outside the match. Downstream code walks the strings with `idx += 1` and indexes into the original polypeptide — truncation shifted those indices by the leading-context offset.
- **P0-2 (crash on empty polypeptides):** `aligner.align("", X)` raises `ValueError`; `pairwise2.align.localds("", X)` returned `[]`. `homoOligomerFixer._alignSeqs` had a `if len(alignments)==0: continue` guard that no longer fired.

Both fixed:
- All 4 alignment files gained a module-level `_padded_local_alignment(seq1, seq2, alignment)` helper.
- `homoOligomerFixer._alignSeqs` gained `if not seq0 or not seq1: return []` at the top.

Byte-for-byte sanity check against legacy pairwise2 on `s1="PPPACDEFGHIKLMNQQQ", s2="ACDEFKLMN"` matched: `'PPPACDEFGHIKLMNQQQ'` / `'---ACDEF---KLMN---'`.

---

## 3. Phase A — env + module loads (complete)

**Cluster Python env**: `protein` (existing conda env at `~/miniforge3/envs/protein`).
**Python version: 3.10.20**.

**Packages installed:**
```bash
conda install -n protein -c conda-forge -c bioconda \
    xgboost gemmi tqdm requests mmtf-python psutil joblib -y
```

Versions on cluster:
- xgboost **3.2.0** — see caveat below
- biopython 1.86, pandas 2.3.3, numpy <2.0, gemmi (≥0.6)

> **xgboost 3.2.0 caveat (still pending):** BIPSPI's `trainAndTest/classifiers/xgBoost.py` was written against xgboost 0.80 (2017). xgboost 3.x has substantial API changes. **Will almost certainly break at training time.** Either pre-empt with `conda install -n protein "xgboost>=1.7,<2.0" -y`, or patch `xgBoost.py` when it errors (per "fix later" scope).

**Module loads** (must repeat in any fresh shell):
```bash
conda activate protein
module load BLAST+/2.14.1-gompi-2023a
module load CD-HIT/4.8.1-GCC-12.2.0
```

## 4. Phase B — al2co + clustalw (complete)

Built/installed via `tools/install_al2co_clustalw.sh`:
- **clustalw 2.1** via bioconda at `~/miniforge3/envs/protein/bin/clustalw`
- **al2co** from [TheApacheCats/al2co GitHub mirror](https://github.com/TheApacheCats/al2co) at `~/tools/al2co/al2co` (with the char[500]→char[1024] buffer patch)

## 5. Phase B.5 — `dependencies.cfg` patcher (complete)

`tools/patch_dependencies_cfg.sh` rewrote `configFiles/cmdTool/dependencies.cfg` with cluster-specific paths:
```
psiBlastBin       /cvmfs/.../BLAST+/2.14.1-gompi-2023a/bin/psiblast
psiBlastDB_path   /home/biostruct01/databases/uniref90/uniref90.fasta
cdHitBin_path     /cvmfs/.../CD-HIT/4.8.1-GCC-12.2.0/bin/cd-hit
clustalW_path     /home/biostruct01/miniforge3/envs/protein/bin/clustalw
al2coBin_path     /home/biostruct01/tools/al2co/al2co
spider2PyScript_path /home/biostruct01/tools/SPIDER2/misc/pred_pssm.py
```

## 6. Phase C — SPIDER2 port (complete)

**Decision taken: option a-i** (numpy-only rewrite). SPIDER2 source: `SPIDER2_local.tgz` from `http://183.36.5.251:8080/sparks_downloads/.../old_versions/SPIDER2_local.tgz` (Zhou lab successor site since Sparks Lab moved). Extracted to `~/tools/SPIDER2/`.

**Key discovery**: `pred_pssm.py` is **already pure numpy** — no Theano needed. The "NN" is matrix multiplies + sigmoid, weights load via `numpy.load()` on bundled `.npz` files.

Three Py2 → Py3 edits applied by `tools/port_spider2.sh`:
1. **3× `print >>fp, X` → `print(X, file=fp)`**
2. **`numpy.load(f)` → `numpy.load(f, allow_pickle=True, encoding='latin1')`** for Py2-pickled `.npz` files

**Validation:** smoke test against bundled `SPIDER2/ex/1a1xA.pssm` produced `1a1xA.spd3` that byte-matches the reference `1a1xA_CHECK.spd3`.

## 7. Phase D — uniref90 BLAST DB (in progress, restarting)

**Status as of 2026-06-03:**
- ✅ uniref90.fasta.gz downloaded (44 GB compressed)
- ✅ gunzip'd to 84 GB FASTA at `~/databases/uniref90/uniref90.fasta`
- ⚠️ `makeblastdb` first attempt started 2026-05-27 18:14, ran for 75 min, was interrupted at 19:28 when terminal closed mid-indexing of volume 19. Result: 19 of ~20 volumes complete, no `.pal` alias file → unusable.
- ⏳ `makeblastdb` restarted 2026-06-03 13:11 inside tmux + srun. Overwrites existing volumes from scratch. ETA ~75 min from restart.

Recovery procedure used:
```bash
tmux new -s uniref90_recovery
srun --partition=cpu --cpus-per-task=4 --mem=16G --time=4:00:00 --pty bash -l
conda activate protein
module load BLAST+/2.14.1-gompi-2023a
cd ~/databases/uniref90/
rm -f uniref90.fasta.19.*               # partial volume
rm -f uniref90.fasta.{pdb,ptf}-lock     # stale locks
rm -f uniref90.fasta.{pdb,ptf,pal}      # partial stubs
cd ~/BIPSPI-Resurrect/BIPSPI
bash tools/download_uniref90_db.sh      # idempotent — skips already-downloaded fasta
# Ctrl-b d to detach. THIS TIME don't close the terminal until detached.
```

When complete, verify with:
```bash
ls ~/databases/uniref90/uniref90.fasta.pal
echo -e ">test\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF" > /tmp/test.fa
psiblast -query /tmp/test.fa -db ~/databases/uniref90/uniref90.fasta -num_iterations 1 -num_threads 4 -out /tmp/test.psiblast
head -30 /tmp/test.psiblast
```

## 8. Runtime fixes (post-port whack-a-mole — three rounds, all applied)

After Phase A-C done, ran `generateBIPSPIModel.py` against bundled example data. Errors surfaced one-by-one, each a small fix:

### Round 1: Biopython 1.80 removals (3 files)
Biopython 1.80 removed `three_to_one` and `one_to_three` from `Bio.PDB.Polypeptide`. Patched with `IUPACData`-based shim functions in:
- `computeFeatures/toolManagerGeneric.py` (uses both)
- `pythonTools/extractSequences.py` (`three_to_one` only)
- `evaluation/workers/showPymolPath.py` (`one_to_three` only)

Shim preserves the original `KeyError`-on-unknown behavior that BIPSPI's `threeLetterAA_to_one` relies on for special-case mappings (MSE→M, PTR→T, etc.). `seq1`/`seq3` would silently swallow these to "X".

### Round 2: Py2/Py3 gzip mode mismatch (7 sites in 5 files)
Py2 `gzip.open(name, "w")` accepted strings; Py3 defaults to binary mode and requires bytes. Patched explicit text mode (`"wt"` / `"rt"`) in:
- `computeFeatures/toolManagerGeneric.py:243` (write)
- `trainAndTest/resultsManager.py:329,367` (writes)
- `computeFeatures/computeFeatsOneComplex.py:107,120` (read + write)
- `utils.py:116,118` (reads — `openForReadingFnameOrGz` helper)
- `pythonTools/myPDBParser.py:37` (read)
- `patchDock/launchPatchDock.py:26` left alone — explicitly `"rb"` for a bytes copy, correct as-is.

### Round 3: gemmi ≥0.6 API change in our adapters (2 files)
gemmi renamed `Model.name` (str) → `Model.num` (int). Patched in:
- `tools/fetch_bipspi_training_set.py`
- `tools/prepare_bipspi_inputs.py`

Both adapters now do `gemmi.Model(str(src_model.num))`.

### Still pending (will surface in next runs)
- `xgboost 0.80` → 3.x API at training time (callbacks, early stopping, etc.) — pre-empt with downgrade or patch on error.
- `from Bio.PDB.Polypeptide import aa1` in `spider2Manager.py` and `Al2coManager.py` — `aa1` constant may also have been removed in 1.80 cleanup. Will likely error when seq feature pipeline runs.

## 9. Canonical training data (NEW — assembled 2026-06-03)

Co-researcher provided `info_HEDt.tab` — BIPSPI's canonical training list with SCOPe-family cluster IDs for K-fold CV grouping:
- 228 uppercase entries (Protein Docking Benchmark 5 — evaluated test complexes)
- 2403 lowercase entries (RCSB augmentation — training-only)
- 5 columns: `pdbId chainIds1 chainIds2 scopes1 scopes2`
- Located at `docs/info_HEDt.tab` in the repo

**Fetcher: `tools/fetch_bipspi_training_set.py`**
- Parses info_HEDt.tab (tolerates header)
- For uppercase entries: copies/symlinks 4-file BIPSPI bundle from `--benchmark5-dir`
- For lowercase entries: parallel RCSB fetch + gemmi chain extraction + `_b → _u` symlinks
- Writes: `pdbs/`, `info_HEDt.noheader.tab` (header stripped for `--scopeFamiliesFname`), `manifest.json`, `failed.txt`
- Idempotent + resumable

**Benchmark 5.5** downloaded from `https://zlab.wenglab.org/benchmark/benchmark5.5.tgz` (45 MB), extracted to `~/benchmark5/benchmark5.5/structures/` (flat dir of `PDBID_{l,r}_{b,u}.pdb` files, BIPSPI's exact naming).

**Full fetch run (2026-06-03):**
- 228 Benchmark 5 entries: 100% (instant, local I/O)
- 2403 RCSB entries: 2402/2403 in 3:13 min (1 failure: `4v4c` — yeast ribosome, too large for legacy PDB format, RCSB only serves it as mmCIF)
- 2610 / 2611 unique complexes prepared (2631 input entries → 2611 unique pdb_ids; 20 duplicate IDs collapsed, matching BIPSPI's one-file-per-pdb convention)
- Output: `~/bipspi_run/canonical/pdbs/` (~10,440 PDB files), `~/bipspi_run/canonical/info_HEDt.noheader.tab`

The single failure (4v4c) is well within noise (0.04%). Not worth a `.cif` fallback for now.

---

## 10. Two paths for the actual training

### Path A (canonical — validate port reproduces published BIPSPI numbers)

Uses the data the canonical fetcher produced. **Run this FIRST** — without confirming the port reproduces published numbers, any score we get on our splits is suspect.

```bash
# After Phase D + bundled-example smoke pass
tmux new -s bipspi_canonical
srun --partition=cpu --cpus-per-task=16 --mem=64G --time=12:00:00 --pty bash -l
conda activate protein
module load BLAST+/2.14.1-gompi-2023a CD-HIT/4.8.1-GCC-12.2.0
cd ~/BIPSPI-Resurrect/BIPSPI
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/canonical/pdbs \
    --scopeFamiliesFname ~/bipspi_run/canonical/info_HEDt.noheader.tab \
    --N_KFOLD 10 \
    --wdir ~/bipspi_run/canonical/wdir \
    --ncpu 16
# Ctrl-b d to detach (this time DON'T close the terminal early)
```
ETA: ~5-8 hours (dominated by 2611 × 2 chains × ~1.5 min psiblast queries against uniref90, parallelised across 16 cores).

### Path B (our splits — apples-to-apples vs 07b RF)

After Path A validates the port. Uses `tools/prepare_bipspi_inputs.py` to convert our existing `~/ESMFold-multimer/data/{splits,structures}/` into BIPSPI format. Single-fold `folds.json` via `--N_KFOLD` instead of `--scopeFamiliesFname`.

```bash
# After Path A validated
python tools/prepare_bipspi_inputs.py \
    --splits-dir ~/ESMFold-multimer/data/splits \
    --structures-dir ~/ESMFold-multimer/data/structures \
    --output-dir ~/bipspi_run/our_splits \
    --evaluate test
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/our_splits/pdbs \
    --N_KFOLD ~/bipspi_run/our_splits/folds.json \
    --wdir ~/bipspi_run/our_splits/wdir \
    --ncpu 16
```

### Open decisions
- **Train+val vs train-only as training set:** adapter defaults to `train+val combined as training, test held out`.
- **Metric harmonisation:** BIPSPI emits `auc_pair`, `prec_50/100/500`, `auc_l/r`, `mcc_l/r`. 07b emits `pair-AUROC/AUPRC`, `prec@L/L2/L5`, `site-AUROC/F1/MCC`. Cleanest comparable: run 07b's evaluator over BIPSPI's per-complex `prefix.tab.res.gz` outputs. Defer until smoke passes.

---

## 11. Reusable inputs from the ESMFold-multimer pipeline (for Path B)

| Pipeline artifact | Used by BIPSPI? | How |
|---|---|---|
| `~/ESMFold-multimer/data/splits/{train,val,test}.json` | Yes | Read directly by `tools/prepare_bipspi_inputs.py` |
| `~/ESMFold-multimer/data/structures/{PDB}_assembly1.cif` | Yes (with conversion) | Adapter extracts chain_a and chain_b as separate PDBs using gemmi |
| `~/ESMFold-multimer/data/labels/*_assembly1.json` | No | BIPSPI computes its own contact maps from the bound PDB |
| `~/ESMFold-multimer/data/cached/esm2.zarr` | No | BIPSPI uses PSI-BLAST PSSMs, not ESM2 embeddings |

---

## 12. Cluster details

- **Login**: `ssh biostruct01@10.205.10.23` (VPN-in-WSL required)
- **Nodes**: `haskell`, `rust` (sibling nodes sharing `/home`)
- **No sbatch** for regular users. Long jobs need `srun --pty bash -l` wrapped in `tmux`.
- **Critical lesson learned 2026-06-03:** `Ctrl-b d` to detach tmux BEFORE closing any terminal. The first Phase D attempt was killed because the terminal was closed without detaching, losing 75 min of makeblastdb progress.

## 13. Files at a glance

```
~/BIPSPI-Resurrect/BIPSPI/                       (cluster)
E:\BIPSPI-Resurrect\BIPSPI\                      (local Windows)
├── STATUS.md                          <- this file
├── bipspi_py3_environ.yml             <- Python 3 conda env spec (reference; we extended `protein` instead)
├── bipspi_plus_environ.yml            <- ORIGINAL Py2.7 env (kept for reference)
├── tools/
│   ├── prepare_bipspi_inputs.py       <- our-splits → BIPSPI pdbsIndir (Path B adapter)
│   ├── fetch_bipspi_training_set.py   <- info_HEDt.tab → BIPSPI pdbsIndir (Path A canonical fetcher)
│   ├── check_cluster_deps.sh          <- haskell dep probe (already run)
│   ├── install_al2co_clustalw.sh      <- Phase B installer (already run, idempotent)
│   ├── patch_dependencies_cfg.sh      <- cfg path patcher (already run, idempotent)
│   ├── download_uniref90_db.sh        <- Phase D downloader (restarting in tmux now)
│   └── port_spider2.sh                <- Phase C SPIDER2 patcher (already run, idempotent)
├── docs/
│   ├── info_HEDt.tab                  <- canonical BIPSPI training list (2631 entries with scopes)
│   ├── scopes_example.tab             <- BIPSPI's 2-complex bundled example
│   └── trainingPDBsExample/           <- BIPSPI's 2-complex bundled PDBs
├── configFiles/cmdTool/
│   └── dependencies.cfg               <- PATCHED with cluster paths (backup at dependencies.cfg.bak.*)
├── Config.py                          <- PORTED
├── generateBIPSPIModel.py             <- (unchanged top-level entry)
├── predictComplexes.py                <- PORTED
├── monitorScreenlog.py                <- PORTED (with one pre-existing parse error left intact)
├── computeFeatures/
│   ├── toolManagerGeneric.py          <- PORTED + three_to_one/one_to_three shim + gzip mode fix
│   ├── computeFeatsOneComplex.py      <- gzip mode fix
│   └── common/{boundUnboundMapper,homoOligomerFixer}.py    <- PORTED + critic P0-1/P0-2 fixes
├── codifyComplexes/
│   ├── codifyProtocols/DataLoaderClass.py     <- PORTED
│   └── codifyProtocols/SeqProtocol.py         <- (unchanged)
├── trainAndTest/
│   ├── trainAndTest.py                        <- PORTED
│   ├── resultsManager.py                      <- PORTED + gzip mode fix
│   └── evaluateResults.py                     <- PORTED
├── pythonTools/
│   ├── alignSequences.py                      <- PORTED + critic P0-1 fix
│   ├── extractSequences.py                    <- + three_to_one shim
│   ├── extractModelsFromPdbFile.py            <- PORTED
│   ├── combinePDBs.py                         <- PORTED
│   └── myPDBParser.py                         <- gzip mode fix
├── utils.py                                   <- gzip mode fix (openForReadingFnameOrGz)
└── evaluation/                                <- 6 files PORTED (.ix/.append idiom)
    └── workers/showPymolPath.py               <- + one_to_three shim

External (on cluster, NOT in this repo):
~/tools/al2co/al2co                            <- al2co binary (Phase B)
~/tools/SPIDER2/                               <- SPIDER2 source + .npz weights (Phase C)
~/tools/SPIDER2/misc/pred_pssm.py              <- PORTED via tools/port_spider2.sh
~/databases/uniref90/uniref90.fasta            <- 84 GB raw FASTA + BLAST DB (Phase D restarting)
~/miniforge3/envs/protein/                     <- conda env (Phase A extended)
~/benchmark5/benchmark5.5/structures/          <- Benchmark 5.5 PDBs (Path A data, 228 complexes × 4 files)
~/bipspi_run/canonical/                        <- canonical Path A inputs ready (2610 complexes)
```

`git diff HEAD` shows every line of the port relative to upstream.
