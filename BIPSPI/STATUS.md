# BIPSPI Python 3 Port — Status & Cluster Handoff

**Upstream:** [rsanchezgarc/BIPSPI](https://github.com/rsanchezgarc/BIPSPI) (Python 2.7)
**Our fork:** [solsylph/BIPSPI-Resurrect](https://github.com/solsylph/BIPSPI-Resurrect)
**Purpose:** sequence-mode resurrection on the ESMFold-multimer mmseqs2-clustered splits, as a non-RF baseline after the RF-on-ESM2 architectural ceiling was confirmed (pair-AUROC ~0.55 on test).
**Last update:** 2026-06-03 (Phase D complete, seq-mode smoke in progress).

> **Repo layout:** GitHub repo has `BIPSPI/` as a subdirectory at its root (local git repo was reinitialised at `E:\BIPSPI-Resurrect\` instead of `E:\BIPSPI-Resurrect\BIPSPI\`). On the cluster: `~/BIPSPI-Resurrect/BIPSPI/`. All commands assume `cd ~/BIPSPI-Resurrect/BIPSPI` first.

---

## Progress at a glance

| Phase | Subject | Status |
|---|---|---|
| Port | Py2.7 → Py3 mechanical port (21 files) | ✅ |
| Critic | Two P0 regressions found + fixed | ✅ |
| A | Extend `protein` conda env + module-load BLAST+/CD-HIT | ✅ |
| B | Build al2co + install clustalw 2.1 (bioconda) | ✅ |
| Cfg patch | Update `dependencies.cfg` with cluster paths | ✅ |
| C | Port SPIDER2 to Py3.10 (numpy-only, option a-i) | ✅ |
| **D** | **Download uniref90 + makeblastdb** | **✅ (2026-06-03)** |
| Runtime fixes | 4 rounds of post-port whack-a-mole (Biopython 1.80, gzip mode, gemmi 0.6, BIPSPI default-arg bug) | ✅ |
| Canonical data | Fetch BIPSPI's published training set (info_HEDt.tab + Benchmark 5 + RCSB) | ✅ (2610/2611 prepared) |
| Smoke seq | `--modelType seq` against `./docs/trainingPDBsExample` | ⏳ **psiblast running 2/6 chains done** |
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

`python -m compileall .` passes for 122 of 123 files. The 1 failure (`monitorScreenlog.py` line 28) is a pre-existing source bug — out of scope.

## 2. Critic/implementor pair-programming round (complete)

Two P0 regressions found and fixed:
- **P0-1**: `PairwiseAligner` local mode returns TRUNCATED aligned strings; `pairwise2.align.localds` returned full-length with `-` padding. Downstream index-walking code was silently producing wrong residue mappings. Fixed via `_padded_local_alignment(seq1, seq2, alignment)` helper in all 4 alignment files.
- **P0-2**: `aligner.align("", X)` raises `ValueError`; `pairwise2.align.localds("", X)` returned `[]`. Fixed via `if not seq0 or not seq1: return []` guard in `homoOligomerFixer._alignSeqs`.

Byte-for-byte sanity check against legacy pairwise2 passed.

---

## 3. Phase A — env + module loads (complete)

**Cluster Python env**: `protein` (existing conda env at `~/miniforge3/envs/protein`).
**Python version: 3.10.20**.

```bash
conda install -n protein -c conda-forge -c bioconda \
    xgboost gemmi tqdm requests mmtf-python psutil joblib -y
```

Versions: xgboost **3.2.0** (caveat below), biopython 1.86, pandas 2.3.3, numpy <2.0, gemmi (≥0.6).

> **xgboost 3.2.0 caveat (still pending):** BIPSPI's `trainAndTest/classifiers/xgBoost.py` was written against xgboost 0.80 (2017). Will likely break at training time. Either pre-empt with `conda install -n protein "xgboost>=1.7,<2.0" -y` or patch when it errors.

**Module loads** (must repeat in every fresh shell):
```bash
conda activate protein
module load BLAST+/2.14.1-gompi-2023a
module load CD-HIT/4.8.1-GCC-12.2.0
```

## 4. Phase B — al2co + clustalw (complete)

Via `tools/install_al2co_clustalw.sh`:
- **clustalw 2.1** via bioconda → `~/miniforge3/envs/protein/bin/clustalw`
- **al2co** from [TheApacheCats/al2co GitHub mirror](https://github.com/TheApacheCats/al2co) → `~/tools/al2co/al2co` (with char[500]→char[1024] buffer patch)

## 5. Phase B.5 — `dependencies.cfg` patcher (complete)

`tools/patch_dependencies_cfg.sh` rewrote cluster paths into `configFiles/cmdTool/dependencies.cfg`.

## 6. Phase C — SPIDER2 port (complete)

**Option a-i** (numpy-only rewrite). Source from `http://183.36.5.251:8080/sparks_downloads/.../old_versions/SPIDER2_local.tgz`. Extracted to `~/tools/SPIDER2/`.

Key: `pred_pssm.py` is **already pure numpy** — no Theano. `tools/port_spider2.sh` applied three Py2→Py3 edits:
1. 3× `print >>fp, X` → `print(X, file=fp)`
2. `numpy.load(f)` → `numpy.load(f, allow_pickle=True, encoding='latin1')` for Py2-pickled `.npz`

Validation: byte-match against bundled reference `1a1xA_CHECK.spd3`. ✅

## 7. Phase D — uniref90 BLAST DB (complete 2026-06-03)

**Final state:** 84 GB raw FASTA + 23 BLAST volumes + `.pal` alias file at `~/databases/uniref90/uniref90.fasta`. 188 million sequences (66 billion residues) queryable. Total install: ~95 GB.

**Took 3 attempts to get right:**

| Attempt | Result | Lesson |
|---|---|---|
| 2026-05-27 | Got to volume 19 of ~20 in 75 min, then **terminal closed without detaching from tmux** → process killed | `Ctrl-b d` to detach BEFORE closing any terminal |
| 2026-06-03 #1 (with -hash_index + -parse_seqids, 16 GB mem) | OOM-killed almost immediately | -parse_seqids on UniRef90 needs >16 GB RAM; either bump mem or drop the flag |
| 2026-06-03 #2 (dropped -parse_seqids, kept -hash_index, 64 GB mem) | `BLAST Database creation error: Duplicate seq_ids are found: GNL|BL_ORD_ID|835307` | -hash_index without -parse_seqids breaks ordinal-ID uniqueness check |
| 2026-06-03 #3 (dropped both flags, 64 GB mem) | ✅ Completed in 35 min | Minimal makeblastdb invocation: `-in / -dbtype / -out` only |

Final `tools/download_uniref90_db.sh` runs the bare minimum:
```bash
makeblastdb -in uniref90.fasta -dbtype prot -out uniref90.fasta
```
`-hash_index` is an optimisation, not a requirement; psiblast queries work fine without it. `-parse_seqids` is for `blastdbcmd` accession lookup, which BIPSPI doesn't use.

**Verification (passed):**
```
PSIBLAST 2.14.1+
Database: uniref90.fasta
           188,848,220 sequences; 66,359,825,357 total letters
Results from round 1
```

## 8. Runtime fixes (post-port whack-a-mole — 4 rounds, all applied)

### Round 1: Biopython 1.80 removals (3 files)
Patched `three_to_one` / `one_to_three` with `IUPACData`-based shim functions in `computeFeatures/toolManagerGeneric.py`, `pythonTools/extractSequences.py`, `evaluation/workers/showPymolPath.py`. Shim preserves the `KeyError`-on-unknown behavior BIPSPI's `threeLetterAA_to_one` relies on for special-case mappings (MSE→M, PTR→T, etc.).

### Round 2: Py2/Py3 gzip mode mismatch (7 sites in 5 files)
Py3 `gzip.open(...)` defaults to binary mode. Patched explicit text mode (`"wt"` / `"rt"`) in `toolManagerGeneric.py`, `resultsManager.py`, `computeFeatsOneComplex.py`, `utils.py`, `myPDBParser.py`. (`patchDock/launchPatchDock.py:26` left as `"rb"` — correct for binary copy.)

### Round 3: gemmi ≥0.6 API change in our adapters (2 files)
gemmi renamed `Model.name` (str) → `Model.num` (int). Patched in `tools/fetch_bipspi_training_set.py` and `tools/prepare_bipspi_inputs.py` to use `gemmi.Model(str(src_model.num))`.

### Round 4 (NEW 2026-06-03): BIPSPI default-argument-evaluation bug
**Pre-existing logic bug** (not Py2→Py3 specific): three function definitions in `generateBIPSPIModel.py` had `methodProtocol=conf.modelType` as a default argument. Python evaluates defaults at function-definition time, BEFORE `parse_args()` runs. So `--modelType seq` was silently ignored and BIPSPI always ran in the configFile default mode (`struct`), which then tried to invoke PSAIA.

Fixed two ways for belt-and-braces:
- `configFile.cfg`: changed default `modelType struct` → `modelType seq` (matches our actual use case)
- `generateBIPSPIModel.py`: changed `computeFeatures`, `codifyStep`, `trainAndTest` to `methodProtocol=None` with `if X is None: X = conf.X` lookups inside the function bodies. CLI overrides now actually work.

### Still pending (will surface in remaining runs)
- `xgboost 0.80 → 3.2 API` at training step (`trainAndTest/classifiers/xgBoost.py`)
- `from Bio.PDB.Polypeptide import aa1` in `spider2Manager.py` and `Al2coManager.py` — `aa1` may also have been removed in 1.80 cleanup

## 9. Canonical training data (assembled 2026-06-03)

Co-researcher provided `info_HEDt.tab` — BIPSPI's canonical training list with SCOPe-family cluster IDs:
- 228 uppercase entries (Protein Docking Benchmark 5 — evaluated test complexes)
- 2403 lowercase entries (RCSB augmentation — training-only)
- 5 columns: `pdbId chainIds1 chainIds2 scopes1 scopes2`
- Stored at `docs/info_HEDt.tab`

**Fetcher: `tools/fetch_bipspi_training_set.py`**
- Parses info_HEDt.tab (tolerates header)
- Uppercase: copies/symlinks 4-file bundle from `--benchmark5-dir`
- Lowercase: parallel RCSB fetch + gemmi chain extraction + `_b → _u` symlinks
- Outputs: `pdbs/`, `info_HEDt.noheader.tab`, `manifest.json`, `failed.txt`
- Idempotent + resumable

**Benchmark 5.5** downloaded from `https://zlab.wenglab.org/benchmark/benchmark5.5.tgz` (45 MB), extracted to `~/benchmark5/benchmark5.5/structures/` (flat dir, BIPSPI's exact naming).

**Run results:**
- 228 Benchmark 5 entries: 100% (instant)
- 2403 RCSB: 2402/2403 in 3:13 min (1 failure: `4v4c` yeast ribosome — too large for legacy PDB format, RCSB only serves as mmCIF; 0.04% loss, ignored)
- 2610/2611 unique complexes prepared (20 duplicate IDs collapsed per BIPSPI's one-file-per-pdb convention)
- Output: `~/bipspi_run/canonical/pdbs/` (~10,440 PDB files)

---

## 10. Two paths for training

### Path A — canonical (run FIRST; validates port reproduces published BIPSPI numbers)

```bash
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
# Ctrl-b d to detach. THIS TIME do not close the terminal early.
```
ETA: ~5-8 hours (dominated by 2611 × 2 chains × ~5 min psiblast queries against uniref90, parallelised across 16 cores).

### Path B — our splits (vs 07b RF; runs after Path A validates the port)

Uses `tools/prepare_bipspi_inputs.py` against `~/ESMFold-multimer/data/{splits,structures}/`. Outputs single-fold `folds.json` for `--N_KFOLD`.

```bash
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

---

## 11. Cluster details + lessons learned

- **Login**: `ssh biostruct01@10.205.10.23` (VPN-in-WSL required). Nodes `haskell` + `rust` share `/home`.
- **No sbatch** for regular users. Long jobs need `srun --pty bash -l` inside `tmux`.
- **`Ctrl-b d` to detach tmux** BEFORE closing any terminal. **Cost us 75 min of makeblastdb work on the first Phase D attempt.**
- **`/tmp` is node-local**, NOT shared. If you `--wdir /tmp/X` on node A you can't see it from node B. **Always use `~/...` paths** for cross-node visibility (the makeblastdb artefacts live in `~/databases/` for this reason). Adapter outputs also use `~/bipspi_run/...`.
- **SSH between compute nodes works** (`ssh haskell` from `rust` prompted for password and connected) — useful for monitoring a long job from a second terminal.
- **Module loads don't persist** across SSH sessions. Re-issue `module load BLAST+/...` and `CD-HIT/...` in every new shell.

## 12. Files at a glance

```
~/BIPSPI-Resurrect/BIPSPI/                       (cluster)
E:\BIPSPI-Resurrect\BIPSPI\                      (local Windows)
├── STATUS.md                          <- this file
├── bipspi_py3_environ.yml             <- Python 3 conda env spec (reference; we extended `protein` instead)
├── bipspi_plus_environ.yml            <- ORIGINAL Py2.7 env (reference)
├── tools/
│   ├── prepare_bipspi_inputs.py       <- our-splits → BIPSPI pdbsIndir (Path B adapter)
│   ├── fetch_bipspi_training_set.py   <- info_HEDt.tab → BIPSPI pdbsIndir (Path A canonical fetcher)
│   ├── check_cluster_deps.sh          <- haskell dep probe (already run)
│   ├── install_al2co_clustalw.sh      <- Phase B installer (already run, idempotent)
│   ├── patch_dependencies_cfg.sh      <- cfg path patcher (already run, idempotent)
│   ├── download_uniref90_db.sh        <- Phase D downloader (run successfully 2026-06-03)
│   └── port_spider2.sh                <- Phase C SPIDER2 patcher (already run, idempotent)
├── docs/
│   ├── info_HEDt.tab                  <- canonical BIPSPI training list (2631 entries with scopes)
│   ├── scopes_example.tab             <- BIPSPI's 2-complex bundled example
│   └── trainingPDBsExample/           <- BIPSPI's 2-complex bundled PDBs
├── configFiles/cmdTool/
│   ├── configFile.cfg                 <- PATCHED: modelType seq (was struct)
│   └── dependencies.cfg               <- PATCHED with cluster paths
├── generateBIPSPIModel.py             <- PATCHED: default-arg eval bug for --modelType CLI overrides
├── Config.py                          <- PORTED
├── predictComplexes.py                <- PORTED
├── monitorScreenlog.py                <- PORTED (with one pre-existing parse error left intact)
├── computeFeatures/
│   ├── toolManagerGeneric.py          <- PORTED + three_to_one/one_to_three shim + gzip mode fix
│   ├── computeFeatsOneComplex.py      <- gzip mode fix
│   └── common/{boundUnboundMapper,homoOligomerFixer}.py    <- PORTED + critic P0-1/P0-2 fixes
├── codifyComplexes/                   <- PORTED
├── trainAndTest/
│   ├── trainAndTest.py                <- PORTED
│   ├── resultsManager.py              <- PORTED + gzip mode fix
│   └── evaluateResults.py             <- PORTED
├── pythonTools/
│   ├── alignSequences.py              <- PORTED + critic P0-1 fix
│   ├── extractSequences.py            <- + three_to_one shim
│   ├── extractModelsFromPdbFile.py    <- PORTED
│   ├── combinePDBs.py                 <- PORTED
│   └── myPDBParser.py                 <- gzip mode fix
├── utils.py                           <- gzip mode fix
└── evaluation/                        <- 6 files PORTED (.ix/.append idiom)
    └── workers/showPymolPath.py       <- + one_to_three shim

External (on cluster, NOT in this repo):
~/tools/al2co/al2co                    <- al2co binary (Phase B)
~/tools/SPIDER2/                       <- SPIDER2 source + .npz weights (Phase C)
~/tools/SPIDER2/misc/pred_pssm.py      <- PORTED via tools/port_spider2.sh
~/databases/uniref90/uniref90.fasta    <- 84 GB raw FASTA + 23 BLAST volumes + .pal (Phase D)
~/miniforge3/envs/protein/             <- conda env (Phase A extended)
~/benchmark5/benchmark5.5/structures/  <- Benchmark 5.5 PDBs (Path A data)
~/bipspi_run/canonical/                <- canonical Path A inputs ready (2610 complexes)
```

`git diff HEAD` shows every line of the port relative to upstream.

---

## 13. Current activity (live)

**As of 2026-06-03 ~15:30**: seq-mode bundled-example smoke test running on haskell, inside tmux+srun.
- Phase D verified ✅ (188M-sequence DB queried successfully)
- Contact maps written ✅
- structComputer correctly skipped ✅ (after Round-4 fix)
- psiblast PSSMs: 2 of 6 chains done (~5 min each at ncpu=2)
- Pending stages: al2co, SPIDER2, codification, xgboost training

Next likely failure: `aa1` import in `spider2Manager.py` / `Al2coManager.py` (Biopython 1.80 removal), OR `xgboost 3.x` API at training step. Each a small targeted patch.

If smoke completes successfully → fire Path A canonical run in fresh tmux+srun, walks itself overnight.
