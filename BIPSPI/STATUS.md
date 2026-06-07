# BIPSPI Python 3 Port — Status & Cluster Handoff

**Upstream:** [rsanchezgarc/BIPSPI](https://github.com/rsanchezgarc/BIPSPI) (Python 2.7)
**Our fork:** [solsylph/BIPSPI-Resurrect](https://github.com/solsylph/BIPSPI-Resurrect)
**Purpose:** sequence-mode resurrection on the ESMFold-multimer mmseqs2-clustered splits, as a non-RF baseline after the RF-on-ESM2 architectural ceiling was confirmed (pair-AUROC ~0.55 on test).
**Last update:** 2026-06-07 (Round 19 applied — corrupt joblib pickle skip).

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
| Runtime fixes | **16 rounds of post-port whack-a-mole** (Rounds 1-10 see below + Rounds 11-16: GroupKFold n_sub<2, md5 encode, zip subscript, pandas numeric_only, xgboost label encoding, mergeSplitFolds fallback) | ✅ |
| Canonical data | Fetch BIPSPI's published training set (info_HEDt.tab + Benchmark 5 + RCSB) | ✅ (2610/2611 prepared) |
| Smoke seq | `--modelType seq` against `./docs/trainingPDBsExample` | ✅ **Complete 2026-06-03. Both seq + seq_2 stages trained and evaluated. 6 complexes, 16 rounds of fixes.** |
| Canonical run | Full 2611-complex BIPSPI seq training | ⏳ **RUNNING** (2026-06-04, restarted with 128 CPUs after 16-CPU run was too slow; ~15 hr ETA for features) |
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

### Round 4 (2026-06-03): BIPSPI default-argument-evaluation bug
**Pre-existing logic bug** (not Py2→Py3 specific): three function definitions in `generateBIPSPIModel.py` had `methodProtocol=conf.modelType` as a default argument. Python evaluates defaults at function-definition time, BEFORE `parse_args()` runs. So `--modelType seq` was silently ignored and BIPSPI always ran in the configFile default mode (`struct`), which then tried to invoke PSAIA.

Fixed two ways for belt-and-braces:
- `configFile.cfg`: changed default `modelType struct` → `modelType seq` (matches our actual use case)
- `generateBIPSPIModel.py`: changed `computeFeatures`, `codifyStep`, `trainAndTest` to `methodProtocol=None` with `if X is None: X = conf.X` lookups inside the function bodies. CLI overrides now actually work.

### Round 5 (2026-06-03): parsePsiBlast None comparison
**Py2 silently treated `None >= int` as False; Py3 raises TypeError.** `parsePsiBlast` initialised `identity, evalue = None, None` then on the first `>` line in the BLAST output tried the validity check before any `Identities = ...` line had been parsed. Fix: `if identity is not None and evalue is not None and ...` guards at both validity-check sites in `al2coWorkers/parsePsiBlast.py`.

### Round 6 (2026-06-03): zip-as-list in seqToolManager
**Py3 `zip()` returns an iterator, not a list.** `list_of_profilesList = zip(*feats)` was later indexed with `list_of_profilesList[featNum][i]` → TypeError. Fix: `list_of_profilesList = list(zip(*feats))` in `seqToolManager.py:55`.

### Round 7 (2026-06-03): Popen bytes-vs-str in subprocess output checks
**Py3 `Popen.communicate()` returns bytes; Py2 returned str.** Four tool-managers had checks like `output[1] != ""` which evaluate True in Py3 because `b"" != ""`. Spuriously triggered "error" path even when subprocess succeeded. Patched with bytes literals (`b""`) + `len() > 0` + `.decode()` where downstream uses str ops:
- `Al2coManager.runCdHit` (cd-hit invocation)
- `Al2coManager.runClustalW` (clustalw invocation)
- `psaiaManager.computeOneFile` (PSAIA — inactive in seq mode, preventive)
- `CCMPredManger`, `PSICOVManager` (inactive in seq mode, preventive)

### Round 8 (2026-06-03): clustalw 2.1 header incompatibility with al2co
**Modern clustalw 2.1 writes `CLUSTAL 2.1 Multiple Sequence Alignments` as its first line; al2co's parser only recognises `CLUSTAL W` (1.83 format).** al2co then treated the title as a sequence name and errored with `Names do not match, was: CLUSTAL, now: InputSeq`. Fix: `Al2coManager.runClustalW` now rewrites the first line of the clustalw output to start with `CLUSTAL W ` before passing to al2co. One-shot post-process, no recompile needed.

### Round 9 (2026-06-03): Biopython model-id assumption in build_peptides
**Biopython's `PPBuilder.build_peptides(structure)` assumes the first model has id == 0.** Some PDBs preserve a non-zero MODEL number → `KeyError: 0`. Fix: `computeContactMap.build_peptides` now bypasses Structure/Model levels entirely by iterating Chains directly and calling `self.ppb.build_peptides(chain, aa_only=False)` on each Chain (level "C"), which has no model-id assumption. Also handles empty Structure → returns `[]` (avoids PEP 479 `StopIteration → RuntimeError` inside joblib generators).

### Round 10 (2026-06-03): empty alignment table in build_correspondence
**`boundUnboundMapper.build_correspondence` called `np.max(aligU2BScores)` on a potentially zero-size array** (when bound or unbound chain list was empty due to malformed PDB). Fix: early return when `aligU2BScores.size == 0`. Makes the run resilient to per-complex data quality issues; downstream treats empty `boundToUnboundDict` as "no correspondence".

### Bundled-data quality issues discovered + fix-up script
Several `*_u.pdb` files in `./docs/trainingPDBsExample/` are **12-byte placeholders** containing just the text `"<prefix>_<l|r>_b.pdb"` instead of being proper symlinks. Biopython's PDBParser sees zero atoms and returns an empty Structure → contact-map and bound/unbound mapping crash downstream.

Confirmed affected: **2c1o, 2v6x**. Likely others in the bundled dir; same pattern may appear in real Protein Docking Benchmark entries fetched from elsewhere.

**Fix-up script** `tools/fix_bundled_data_stubs.sh` finds every `*_u.pdb` under 100 bytes in a given directory and replaces it with a proper `ln -s <prefix>_<l|r>_b.pdb <prefix>_<l|r>_u.pdb` symlink, matching BIPSPI's documented convention (`docs/repo_help.md`). Idempotent.

```bash
# preemptive sweep on the bundled example dir:
bash tools/fix_bundled_data_stubs.sh
# or on any directory:
bash tools/fix_bundled_data_stubs.sh /path/to/pdbsIndir
# script will also print the cache-cleanup commands you may want to run
# afterwards if you're re-using an existing wdir.
```

**Note:** these are runtime modifications to the local filesystem on the cluster — they're **not** committed to git (would diverge from BIPSPI upstream). Re-run the fix-up script after any fresh clone or after pulling new bundled data.

Also, **`computeFeatures/common/boundUnboundMapper.py:build_correspondence`** has an empty-array guard (Round 10 above), so even without the fix-up script, malformed complexes now skip gracefully instead of crashing the whole batch. The script is the proactive solution; the guard is the safety net.

### Round 11 (2026-06-03): GroupKFold n_sub < 2 in crossValidationSplitter
**`GroupKFold` requires `n_splits >= 2`; with 6 example complexes and N_KFOLD=2, each outer training split falls into 1 scope group.** Two fixes: (a) cap `n_sub = min(N_SUB_FOLDS, len(set(current_groups)))` before calling `GroupKFold`; (b) when `n_sub < 2`, fall back to a single sub-fold using all training indices (degenerate case — canonical run always has `n_groups >= 3`). Also fixed `e.message` → `str(e)` (Py2→Py3) in the except handler.

### Round 12 (2026-06-03): hashlib.md5 requires bytes (Py2→Py3)
**`hashlib.md5(str)` worked in Py2 (str was bytes); Py3 requires encoded bytes.** Fixed `hashlib.md5("...".encode())` in `trainAndTest/processOneFold.py:156`.

### Round 13 (2026-06-03): xgboost 3.x label encoding
**xgboost 3.x `binary:logistic` requires labels `{0, 1}`; BIPSPI uses `{-1, 1}`.** Fixed by remapping: `trainLabels = (np.asarray(trainLabels) > 0).astype(np.float32)` before `modelo.fit()` in `trainAndTest/classifiers/xgBoost.py`. `predict_proba` output is unchanged.

### Round 14 (2026-06-03): zip subscript in processOneFold (Py2→Py3)
**`zip(*list)[1]` fails in Py3 — zip returns an iterator.** Fixed `list(zip(*resultsForEvaluation_list))[1]` in `trainAndTest/processOneFold.py:201`.

### Round 15 (2026-06-03): pandas 2.x mean() on mixed-type DataFrame
**`df.mean(axis=0)` in pandas 2.x raises `TypeError` on non-numeric columns (the prefix-name column).** Fixed `summary.mean(axis=0, numeric_only=True)` in `trainAndTest/trainAndTest.py:225`.

### Round 16 (2026-06-03): mergeSplitFolds stage-2 orthogonality fallback
**With 6 example complexes and N_KFOLD=2, every stage-1 prediction was made by a model trained on exactly the stage-2 test set — orthogonality is mathematically impossible.** Added a fallback in `mergeSplitFolds`: when `trainPrefixes_idx` is empty, warn and use all stage-1 predictions for training complexes regardless of overlap. Raises only if still empty after fallback. Canonical run (2610 complexes, N_KFOLD=10) never hits this path.

### Round 17 (2026-06-04): pd.concat([]) guard in DataLoaderClass + codification skip
**`pd.concat([])` raises `ValueError: No objects to concatenate` in pandas 2.x** when a complex has no loadable feature files (e.g., al2co emits "No alignments read" for very short chains). Two fixes: (a) guard in `codifyComplexes/codifyProtocols/DataLoaderClass.py` — raise `ValueError` before concat if `resultDF_list` is empty; (b) wrap `launchCodifyOneComplex` body in try-except to skip and warn instead of crashing the whole batch.

### Round 18 (2026-06-04): getScopeGroups skips short lines
**`info_reduced.tab` had two lines with an empty `chainIds2` field (double tab → only 4 tokens on split).** Caused `ValueError: not enough values to unpack (expected 5, got 4)` in `loadFromTable`. Fix: skip lines with `len(lineArray) < 5` with a warning in `trainAndTest/getScopeGroups.py`. (The two offending lines were duplicate entries for 1OYV and 1QFW; valid entries exist on adjacent lines.)

### Round 19 (2026-06-07): corrupt joblib pickle skip in predictOnePrefix
**A joblib pickle for a test complex was truncated** (written incompletely when a previous run was killed mid-write), causing `ValueError: EOF: reading array data, expected 262144 bytes got 52928` in `getDataForTestFromPrefix`. Fix: `predictOnePrefix` now delegates to `_predictOnePrefixInner` and wraps it in try-except, returning `None` on failure with a warning. Caller (`trainAndTestOneFold`) filters `None` values from `resultsForEvaluation_list` before processing. Consistent with Round 17 pattern in `launchCodifyOneComplex`. Corrupt complexes are skipped and evaluation continues.

### Confirmed NOT issues
- `from Bio.PDB.Polypeptide import aa1` still works in Biopython 1.86.
- `nthread` param in xgboost 3.x — deprecated but accepted without error.

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

**As of 2026-06-03**: **Smoke test COMPLETE.** 16 rounds of runtime fixes total (Rounds 1-10 feature pipeline; Rounds 11-16 training pipeline). Both seq and seq_2 stages ran end-to-end on the 6-complex bundled example. Results produced (auc_pair ~0.61 stage-1, ~0.59 stage-2 — meaningless for 6 toy complexes, proves the pipeline runs). Models saved at `/tmp/test_bipspi_seq_v2/modelsComputed/`.

**As of 2026-06-04**: **Path A canonical run LIVE** on haskell in `tmux attach -t bipspi_canonical`.

First attempt ran with 16 CPUs (8 parallel psiblasts) — projected ~13 days, too slow. Killed after ~12 hours (~234 chains cached). Restarted with 128 CPUs (64 parallel psiblasts) and 7-day time limit.

Current srun:
```bash
srun --partition=cpu --cpus-per-task=128 --mem=128G --time=7-00:00:00 --pty bash -l
```

**ETA**: ~15 hours for feature computation (128 CPUs = 64 parallel psiblasts, ~5800 remaining chains). Training after that. All previously computed chains (~234) are cached and will be skipped.

**To check progress**:
```bash
DONE=$(ls ~/bipspi_run/canonical/wdir/computedFeatures/common/contactMaps/*.cMap.tab.gz 2>/dev/null | wc -l)
echo "$DONE / 2610 ($(( DONE * 100 / 2610 ))%)"
```

When you see `All features computed for: 2610` the slow part is done.

**After canonical run**: check auc_pair on the Benchmark 5 uppercase complexes against BIPSPI's published ~0.72, then fire Path B.

**Path B** (our mmseqs2 splits vs 07b RF — runs after canonical validates the port):
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

### Quick handoff for a fresh chat
1. Open Claude Code with CWD = `E:\BIPSPI-Resurrect\BIPSPI\`
2. Prompt: "Read STATUS.md. Canonical run in progress / done. Last output: [paste]."
3. STATUS.md is the single source of truth. No other files required.
