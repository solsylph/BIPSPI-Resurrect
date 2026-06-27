# BIPSPI Python 3 Port — Status & Cluster Handoff

**Upstream:** [rsanchezgarc/BIPSPI](https://github.com/rsanchezgarc/BIPSPI) (Python 2.7)
**Our fork:** [solsylph/BIPSPI-Resurrect](https://github.com/solsylph/BIPSPI-Resurrect)
**Purpose:** sequence-mode resurrection on the ESMFold-multimer mmseqs2-clustered splits, as a non-RF baseline after the RF-on-ESM2 architectural ceiling was confirmed (pair-AUROC ~0.55 on test).
**Last update:** 2026-06-26 (Path B codification UNBLOCKED — four chained blockers fixed: `@`-prefix collision, ESM2 resId/`setCurrentSeq`, broken `_u.pdb` symlinks, ESM2 `pdbId` lookup. `dbg_codify` → SUCCESS, ESM2 misses back to 13/0.17%. Full `generateBIPSPIModel.py` run is the next step — see §15).

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
| ESM2 integration | Add ESM2 per-token embeddings as BIPSPI features | ✅ implemented (2026-06-10) |
| ESM2 feature compute | Run feature step on our splits (Path B) | ✅ **clean 2026-06-26** (97.6% retained after SEQRES-hash fix; see §15) |
| ESM2 codification | One-complex codify validates on Path-B features | ✅ **2026-06-26** (`dbg_codify 3wg7jl` → SUCCESS after 4 chained fixes; resId join correct; §15) |
| Re-baseline | Our splits + 07b-comparable evaluation | ⏳ **ready to launch 2026-06-26** — full `generateBIPSPIModel.py` run (`--ncpu 16 --mem=192G --time=12:00:00`); pair-AUROC table vs 07b RF ~0.55 is the deliverable (§15) |

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

**As of 2026-06-03**: **Smoke test COMPLETE.** 16 rounds of runtime fixes (Rounds 1–16). Both seq + seq_2 stages ran end-to-end on the 6-complex bundled example.

**As of 2026-06-08**: **Reduced canonical run COMPLETE.** Switched from the full 2610-complex run (too slow at 16 CPUs) to a reduced set (`info_reduced.tab`, ~1647 complexes with pre-cached features). Three additional runtime fix rounds:
- Round 17: `pd.concat([])` guard in `DataLoaderClass.py` + codification skip on error
- Round 18: `getScopeGroups.py` skips short/malformed lines in scope table
- Round 19: corrupt joblib pickle skip in `processOneFold.py`

OOM issues during `seq_2` codification resolved by requesting 256G RAM from SLURM (node has 503G total). Final run used `--ncpu 32 --mem=256G --time=8:00:00`. The model was trained and saved at `~/bipspi_run/canonical/wdir/modelsComputed/model.seq_2` moments before the 8-hour time limit hit.

**Results on 27 Benchmark 5 uppercase complexes** (in the reduced set — not the full 228-complex published benchmark):

| complex | auc_pair | … |
|---|---|---|
| mean (27 complexes) | **0.8485** | see stdout table from the run |

Note: 0.85 is not directly comparable to BIPSPI's published ~0.72 (different complex set, only 27 of 228). The pipeline is working correctly.

**As of 2026-06-10**: **ESM2 integration COMPLETE.** Added `Esm2Manager.py` + config flags to replace classic PSSM/al2co/SPIDER2 features with pre-computed ESM2-650M per-token embeddings. See §14 for full details.

**Next step: Path B** — run BIPSPI with ESM2 features on our mmseqs2-clustered splits for apples-to-apples comparison with 07b RF (pair-AUROC ~0.55).

**Pre-flight before Path B**:
```bash
# On cluster (protein env):
python -c "import zarr; s=zarr.open('/home/biostruct01/ESMFold-multimer/code/esmfold_multimer/data/cached/esm2.zarr','r'); print(len(list(s.keys())), 'entries')"
```

**Path B run** (see §14 for full commands):
```bash
# Enable ESM2 mode in config, then:
python tools/prepare_bipspi_inputs.py \
    --splits-dir ~/ESMFold-multimer/code/esmfold_multimer/data/splits \
    --structures-dir ~/ESMFold-multimer/code/esmfold_multimer/data/structures \
    --output-dir ~/bipspi_run/esm2_splits --evaluate test
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/esm2_splits/pdbs \
    --N_KFOLD ~/bipspi_run/esm2_splits/folds.json \
    --wdir ~/bipspi_run/esm2_splits/wdir \
    --ncpu 32
```

---

## 14. ESM2 Integration (2026-06-10)

### Goal
Replace BIPSPI's classic sequential features (PSSM/al2co/SPIDER2) with pre-computed ESM2-650M per-token embeddings from the ESMFold-multimer pipeline. Allows apples-to-apples comparison with the 07b RF baseline: same proteins, same embeddings, different model architecture (XGBoost pair-coding vs Random Forest).

### Files changed

| File | Change |
|---|---|
| `computeFeatures/seqStep/seqToolManagers/conservationTools/Esm2Manager.py` | **new** — reads `per_tok_embedding [L, 1280]` from zarr, writes per-chain `.esm2.tab.gz` in BIPSPI's standard format |
| `computeFeatures/seqStep/seqFeatsComputer.py` | ESM2 branch + `skipClassicFeatures` gate; `copySameSeq` updated |
| `codifyComplexes/codifyProtocols/SeqProtocol.py` | `FEATURES_TO_INCLUDE_CHAIN_ESM2` + `FEATURES_TO_INCLUDE_PAIR_ESM2` constants |
| `codifyComplexes/codifyOneComplex.py` | Passes ESM2 feature list to `SeqProtocol` when `useEsm2 + skipClassicFeatures` |
| `configFiles/cmdTool/dependencies.cfg` | `useEsm2 False`, `skipClassicFeatures False`, `esm2ZarrPath_path ~/ESMFold-multimer/...` |

### How the zarr is structured
`/home/biostruct01/ESMFold-multimer/code/esmfold_multimer/data/cached/esm2.zarr` — written by `06_cache_esm2.py`:
- Group key: `sha256(sequence)` (SHA256 of the canonical RCSB sequence string)
- Arrays per group: `per_tok_embedding [L, 1280]` float32, `mean_embedding [1280]` float32

### Activation
Enable ESM2-only mode by editing `configFiles/cmdTool/dependencies.cfg`:
```
useEsm2 True
skipClassicFeatures True
```
Or set on cluster before running:
```bash
sed -i 's/^useEsm2 False/useEsm2 True/' configFiles/cmdTool/dependencies.cfg
sed -i 's/^skipClassicFeatures False/skipClassicFeatures True/' configFiles/cmdTool/dependencies.cfg
```

### Path B run command (ESM2 features, our mmseqs2 splits)
```bash
# 1. Prepare our splits as BIPSPI inputs
python tools/prepare_bipspi_inputs.py \
    --splits-dir ~/ESMFold-multimer/code/esmfold_multimer/data/splits \
    --structures-dir ~/ESMFold-multimer/code/esmfold_multimer/data/structures \
    --output-dir ~/bipspi_run/esm2_splits \
    --evaluate test

# 2. Enable ESM2 mode
sed -i 's/^useEsm2 False/useEsm2 True/' configFiles/cmdTool/dependencies.cfg
sed -i 's/^skipClassicFeatures False/skipClassicFeatures True/' configFiles/cmdTool/dependencies.cfg

# 3. Run (no psiblast needed — features are loaded from zarr)
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/esm2_splits/pdbs \
    --N_KFOLD ~/bipspi_run/esm2_splits/folds.json \
    --wdir ~/bipspi_run/esm2_splits/wdir \
    --ncpu 32
```

### Known risk: sequence hash mismatch → CONFIRMED + FIXED (2026-06-26, see §15)
The zarr was built from RCSB **SEQRES** sequences (from `candidates.json`); BIPSPI extracts sequences from PDB **ATOM** records (resolved residues only). The "both use the same CIF so hashes match" assumption was **WRONG** — 1,857/7,609 complexes (~25%) missed on the first real run. Two parts to the fix: (1) the zarr key is `sha256(seq)[:16]` not the full digest (`Esm2Manager` was computing 64 chars — commit `16f7a92`); (2) a **SEQRES-candidate fallback** in `Esm2Manager` that, on an ATOM-hash miss, looks up the candidate SEQRES sequence for the same pdbId, fetches its `[L,1280]` embedding, and global-aligns ATOM→SEQRES to slice one embedding row per resolved residue (commit `d03e479`). This dropped misses from 1,857 → 13 (0.17%).

### Expected outcome
| Model | Features | Arch | Expected pair-AUROC |
|---|---|---|---|
| 07b RF | ESM2 1280d | Random Forest | ~0.55 (measured) |
| BIPSPI + ESM2 | ESM2 1280d | XGBoost + pair-coding | TBD |
| BIPSPI canonical | PSSM + al2co + SPIDER2 | XGBoost | ~0.72 (published) |

If BIPSPI+ESM2 > RF+ESM2: the XGBoost pair-coding architecture matters more than the features.
If BIPSPI+ESM2 ≈ RF+ESM2: the architectural ceiling is in the features themselves.

### Pre-flight checks before Path B run
```bash
# zarr accessible?
python -c "import zarr; s=zarr.open('~/ESMFold-multimer/code/esmfold_multimer/data/cached/esm2.zarr','r'); print(len(list(s.keys())), 'entries')"
# zarr installed?
python -c "import zarr; print(zarr.__version__)"
# splits ready?
ls ~/bipspi_run/esm2_splits/pdbs/ | head -5
```

### Quick handoff for a fresh chat
1. Open Claude Code with CWD = `E:\BIPSPI-Resurrect\BIPSPI\`
2. Prompt: "Read STATUS.md. Canonical run in progress / done. Last output: [paste]."
3. STATUS.md is the single source of truth. No other files required.

---

## 15. Path B run attempt + blockers (2026-06-26)

First real Path-B run on the cluster (`haskell` node, `~/bipspi_run/esm2_splits/`). Feature computation now runs **clean**; codification is **blocked** on an input-naming collision. Four commits landed this session; one fix is staged locally and not yet pushed.

### Commits landed (pushed)
| Commit | Round | Fix |
|---|---|---|
| `16f7a92` | — | `Esm2Manager._seqHash` truncates to `[:16]` to match the zarr key (`06_cache_esm2.py:243` uses `sha256(seq).hexdigest()[:16]`). Was computing full 64-char digest → 100% miss. |
| `185bc5b` | 20 | `BadNumberOfResidues` made picklable (stores `nResidues`/`partnerId` + `__reduce__`). Worker exceptions were unpicklable in the parent (`__init__` needs 2 args, pickle supplied 1 msg string) → killed the multiprocessing `_handle_results` thread → risk of hang at join. |
| `50c9146` | 21 | `launchComputeFeaturesOneComplex` wrapped in try-except → skip bad complex with warning instead of `_raise_error_fast` aborting the whole joblib batch. Same skip-don't-crash pattern as Rounds 17/19 (codify/predict steps); this was the last unguarded step (feature compute). |
| `d03e479` | — | `Esm2Manager` SEQRES-candidate fallback with align+slice (see §14 risk note). Dropped ESM2 misses 1,857 → 13. Adds `esm2CandidatesJson_path` to `dependencies.cfg` + threads it through `seqFeatsComputer.py`. |

### Fixes this session (all verified working on cluster; doc-only edits may be unpushed — user pushes, see memory `git-push-is-users-job`)
Five chained code fixes were needed to get one complex to codify cleanly (each unmasked the next):
- **`SeqProtocol.py`: `FEATURES_TO_INCLUDE_PAIR_ESM2 = None`** (was `[]`). Empty list passed the `if not self.pairfeatsToInclude is None:` guard in `AbstractProtocol.applyProtocol:87` then indexed `self.pairfeatsToInclude[0]` at `addPairFeatures:224` → `IndexError` on every complex → `0 train complexes loaded`. `None` skips `addPairFeatures`.
- **`tools/prepare_bipspi_inputs.py`: `make_prefix` now `@`-free** (`f"{pdb}{chains}"`) — BLOCKER 1 (`@`-collision).
- **`tools/rename_strip_at_prefixes.sh`** (new): rename existing run's artifacts in place (kept ~1 hr ESM2 compute). Step **1b** also re-points `_u.pdb` symlinks — BLOCKER 3.
- **`Esm2Manager.py` (two fixes):** (a) `seqStructMap.setCurrentSeq(...)` before the residue loop — BLOCKER 2 (resId `0?` → real PDB numbers); (b) `pdbId = prefix[:4].lower()` (was `prefix.split("@")[0]`) — BLOCKER 4 (SEQRES-fallback lookup).
- **`dbg_codify.py`** (new): committable single-complex codify smoke test (replaces the lost ad-hoc cluster version).
- **`dbg_compute_feats.py`** (new): runs ONLY the feature-compute step (auto-removes the `allFeaturesComputed.txt` marker; regenerates in place, no codify) so features can be re-validated before the full non-cached codify.
- **`STATUS.md`**: this §15 update (BLOCKERs 1–4, progress checklist, full-run command).

Outcome: `dbg_codify 3wg7jl` → **SUCCESS**; ESM2 misses 13/0.17%; ready for the full run (see Progress checklist at end of §15).

### ⛔ BLOCKER 2 (RESOLVED in code): ESM2 resId bug → "dataset is empty"
After the `@`-rename, codify reached `AbstractProtocol.applyProtocol` but every complex failed `assert allPairsCodified.shape[0]>1, "Error, <prefix> dataset is empty"` (`AbstractProtocol.py:94`). Root cause: `combinePairwiseAndSingleChainFeats` inner-joins the contact map against the single-chain feature files on `(chainId, resId, resName)`. The ESM2 files had **`resId` = `0?`, `1?`, `2?`...** (the `str(seqIx)+"?"` fallback at `Esm2Manager.py:201`) instead of real PDB residue numbers (`1`, `2`, `3` in the cMap), so the join matched 0 rows. (`resName` was fine — 1-letter in both.)

The fallback fired because `seqStructMap.seqToStructIndex()` returned `None` for **every** residue: at `seqStructMapper.py:248` it does `self.seqToRefSeq[(chainType,chainId)]`, which raises `KeyError` (caught → `None`) unless `setCurrentSeq()` has registered that chain. **Every** classic single-chain manager (PsiBlast/Spider2/Al2co/HHblits/windowSeq) calls `seqStructMap.setCurrentSeq(seqStr, chainType, chainId)` right after `getSeq`; the ESM2 manager was the only one that didn't. Fix = add that one call (matches the proven pattern; does not touch the shared `seqToStructIndex`). Systematic — confirmed by 6/6 sampled prefixes failing identically.

**Regeneration required:** the on-disk `esm2/*.esm2.tab.gz` have the wrong resIds baked in, and the manager's `os.path.isfile` short-circuit would keep them. **GOTCHA:** `computeFeaturesAllPdbsOneDir` short-circuits the *entire* step if `<computedFeatsRootDir>/allFeaturesComputed.txt` exists (prints "All features computed for: N" and returns) — deleting the esm2 files alone is NOT enough, you must also delete that marker. `dbg_compute_feats.py` now removes the marker automatically. Procedure (cluster, `protein` env, after `git pull`):
```bash
rm ~/bipspi_run/esm2_splits/wdir/computedFeatures/seqStep/esm2/*.esm2.tab.gz
rm -f ~/bipspi_run/esm2_splits/wdir/computedFeatures/allFeaturesComputed.txt   # or rely on dbg_compute_feats.py
cd ~/BIPSPI-Resurrect/BIPSPI
PYTHONPATH=. python dbg_compute_feats.py --modelType seq \
    --pdbsIndir ~/bipspi_run/esm2_splits/pdbs \
    --N_KFOLD  ~/bipspi_run/esm2_splits/folds.json \
    --wdir     ~/bipspi_run/esm2_splits/wdir --ncpu 16   # only esm2 recomputes; rest cached
# re-validate (now expect real resIds + SUCCESS):
zcat ~/bipspi_run/esm2_splits/wdir/computedFeatures/seqStep/esm2/3wg7jl_l_*.esm2.tab.gz | cut -f1-3 | head
PYTHONPATH=. python dbg_codify.py --wdir ~/bipspi_run/esm2_splits/wdir 3wg7jl
```
Only after `SUCCESS` (and the resId column showing `1,2,3...` not `0?,1?...`) launch the full `generateBIPSPIModel.py` run.

### ⛔ BLOCKER 3 (RESOLVED): rename broke `_u.pdb` symlinks
The regeneration run skipped **every** complex with `FileNotFoundError: .../pdbs/<prefix>_l_u.pdb`. Cause: `prepare_bipspi_inputs.py` creates `_u.pdb` as **relative symlinks** to `_b.pdb` (`os.symlink(src.name, dst)`). `mv` renames a symlink file but NOT its stored target string, so after the `@`-strip rename, `4lvhbc_l_u.pdb` still pointed at the old `4lvh@bc_l_b.pdb` → dangling. Fix = re-point the links (`ln -sfn <prefix>_l_b.pdb <prefix>_l_u.pdb`). `tools/rename_strip_at_prefixes.sh` now does this automatically (step "1b"); for the already-renamed run a one-off loop over `*_u.pdb` repointed them. After repointing, regeneration sees real pdbs and recomputes esm2.

### ⛔ BLOCKER 4 (RESOLVED): `@`-removal broke ESM2 SEQRES-fallback pdbId lookup
After repointing, regeneration ran but ~1844 complexes failed with `ESM2 embedding not found ... no candidate SEQRES sequence for pdbId 4lmsab aligned` — note the **wrong pdbId `4lmsab`** (should be `4lms`). `Esm2Manager._lookupEmbedding` extracted the pdbId for the SEQRES-candidate fallback via `prefix.split("@")[0]`, which **relied on the `@`** to separate the 4-char PDB code from the chain letters. With `@` gone, `"4lmsab".split("@")[0]` returns `"4lmsab"`, not a key in `candidates.json` (keyed by 4-char `4lms`), so the fallback found nothing and every ATOM≠SEQRES complex (the ones needing the fallback, incl. `3wg7jl`) failed; fast-path exact-hash hits still wrote esm2 files. Fix = `pdbId = prefix[:4].lower()` (correct for both naming schemes). Regenerate again: fast-path esm2 files are kept (correct resIds from the setCurrentSeq fix), only the previously-failed fallback complexes recompute.

### Feature-compute run results (clean)
- Splits prepared: **7,359 train / 250 test** = 7,609 complexes (`prepare_bipspi_inputs.py` on `~/ESMFold-multimer/code/esmfold_multimer/data/{splits,structures}`).
- zarr: **12,951 entries**, arrays `per_tok_embedding [L,1280]` + `mean_embedding [1280]`, keyed on `sha256(SEQRES)[:16]`.
- After the SEQRES fallback: **13** `ESM2 embedding not found` (0.17%), **186** total skips (2.4%): 159 `BadNumberOfResidues` (legit too-short/long chains), 11 `IndexError`, 3 `EOFError` (truncated cache from a killed run), 1 `ValueError`. **97.6% retained** — comparable coverage to 07b RF.
- Probe confirmed all 1,764 skipped pdbIds were present in `candidates.json` with chain-seq hash in the zarr (1764/1764) — i.e. no proteins are genuinely missing, purely a SEQRES-vs-ATOM key problem.

### ⛔ BLOCKER 1 (RESOLVED): `@`-prefix codification collision
After the `None` fix, codification still fails because **`prepare_bipspi_inputs.py` names complexes with chain-pair letters after `@`** (e.g. `4lvh@bc`, `4lac@ac`, `4lac@bc`), but BIPSPI uses `@` as a **"same-complex sampling variant" tag** and strips everything after it in **~15 places** (`grep -n 'split("@")\[0\]'`): notably `AbstractProtocol:79` (`raw_prefix`), `DataLoaderClass:69-70` (file lookup), `crossValidationSplitter:52,59` (CV fold grouping), `trainAndTest:316,342` + `processOneFold` (result averaging), `ComplexCodified:112` (sampling). Effects:
- File lookup strips `4lvh@bc` → `4lvh`, but feature files are named `4lvh@bc_l_B_u_.esm2.tab.gz` → mismatch.
- `4lac@ac` and `4lac@bc` both collapse to `4lac` → `more than 1 Contact map for 4lac` error, and would be wrongly averaged as one complex in eval.

**Fix DECIDED + STAGED (2026-06-26, "fix script + rename artifacts"):** make each complex's base prefix unique and `@`-free (`4lvh@bc` → `4lvhbc`). The 4-char PDB-ID boundary guarantees no collisions (`4lac@ac`/`4lac@bc` → distinct `4lacac`/`4lacbc`). Two coordinated changes, both committable so the cluster gets them via `git pull`:
1. **`tools/prepare_bipspi_inputs.py:make_prefix`** patched: `f"{pdb}@{chains}"` → `f"{pdb}{chains}"` (the `@` separator is gone). The script is now correct for all future prepare runs — no recurrence.
2. **`tools/rename_strip_at_prefixes.sh`** (new) renames the *existing* run's artifacts in place to match, so the ~1 hr of ESM2 features is kept rather than recomputed. The transform is purely "delete the `@` char" — `@` appears nowhere in these names/contents except as the prefix separator, so stripping it lands exactly on the new scheme. Dry-run by default; `--apply` to act; refuses to clobber. It renames files+dirs under `pdbs/` and `wdir/computedFeatures/` (`find -depth`), strips `@` from `folds.json`/`manifest.json` (with `.bak`), and deletes stale `*.train.pkl.gz`/`*.predict.pkl.gz`.

Validation driver **`dbg_codify.py`** (new, committable — the previous one was ad-hoc on the cluster and lost): `PYTHONPATH=. python dbg_codify.py --wdir ~/bipspi_run/esm2_splits/wdir 4lvhbc` → expects `SUCCESS` before any full relaunch (codify does NOT cache, so a failed full run wastes the whole encode).

### Co-researcher guidance (2026-06-26) — confirms the integration approach
BIPSPI's author (rsanchezgarc) provided reference links. They **validate our manager approach** and recommend one simplification:
- Build the ESM manager modeled on `PsiBlastManager.py` (see `seqToolManager.py:30` base class, `PsiBlastManager.py:62,240`, `seqStructMapper.py` for seq→struct mapping). ✅ done as `Esm2Manager.py`.
- Add import + instantiate in `seqFeatsComputer.py`. ✅ done.
- **In `SeqProtocol.py`: comment out the classic-feature line `FEATURES_TO_INCLUDE_CHAIN` L15 (`predAsaAndSS`/SPIDER2) and replace with the ESM features directly** — i.e. edit the `FEATURES_TO_INCLUDE_CHAIN` list in place, rather than our parallel `FEATURES_TO_INCLUDE_CHAIN_ESM2` constant + `useEsm2/skipClassicFeatures` branching in `codifyOneComplex.py`. Our branching approach is functionally equivalent and already works; the next chat can keep it or simplify to the in-place edit per the author's pattern. Either way the `@`-collision blocker is orthogonal and must be fixed first.

### Cluster resource notes (this session)
- Node `haskell`: 512 CPU / ~503 G RAM. Node was ~190 G free, 428 CPU idle. QOS caps: `normal` = 128 CPU/user, `student` = 64 CPU/user. Partition `cpu` default 4 h, max 7 d, `OverSubscribe=FORCE:2` (CPU shareable, memory exclusive).
- **`squeue`/`sacct`/`sinfo` live only on the login node `rust`, NOT on compute nodes** — diagnose job state from the run log + `top`/`free` on the compute node.
- **ESM2 mode is memory-bound, not CPU-bound** (psiblast — the old CPU hog — is skipped; zarr reads are ms). `--ncpu` drives both parallelism AND peak codify RAM. Right-sized ask: `--cpus-per-task=16 --mem=64G`, but **codification OOM'd at 64 G near the end of stage-1** (~7,178 complexes). Relaunch needs **`--mem=192G`**. Set `--ncpu` == `--cpus-per-task` (never exceed allocation).
- **Codification does NOT cache per-complex to disk** — a restart re-encodes from scratch (~same time), it does not resume. Feature computation DOES cache (`os.path.isfile` skip), so restarts there are cheap.
- Terminal paste tip: long one-line commands wrap and break at the filename in this SSH/tmux setup. `cd` into the dir and use short relative paths, or use heredocs.

### Progress checklist (2026-06-26)
**DONE on cluster** (`~/bipspi_run/esm2_splits/`, env `protein`):
1. ✅ Pushed the fix batch; `git pull` on cluster.
2. ✅ `@`-rename applied (`tools/rename_strip_at_prefixes.sh ... --apply`): 93,999 files renamed, 0 `@` left, `_u.pdb` symlinks repointed.
3. ✅ ESM2 features regenerated (`dbg_compute_feats.py`): resId column now `1,2,3...` (not `0?,1?...`); ESM2 misses back to **13** (0.17%).
4. ✅ One-complex codify validated: `dbg_codify.py --wdir .../wdir 3wg7jl` → **SUCCESS** (codified in 0.25s).

**NEXT — the only remaining step: full Path-B run.** In `tmux`, on a compute node, detach with `Ctrl-b d` before closing any terminal:
```bash
tmux new -s bipspi_pathB
srun --partition=cpu --cpus-per-task=16 --mem=192G --time=12:00:00 --pty bash -l
conda activate protein
cd ~/BIPSPI-Resurrect/BIPSPI
# pre-flight: ESM2 mode on + feature step is a cached no-op
grep -E '^(useEsm2|skipClassicFeatures)' configFiles/cmdTool/dependencies.cfg   # both True
ls ~/bipspi_run/esm2_splits/wdir/computedFeatures/allFeaturesComputed.txt        # exists -> compute skips
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/esm2_splits/pdbs \
    --N_KFOLD  ~/bipspi_run/esm2_splits/folds.json \
    --wdir     ~/bipspi_run/esm2_splits/wdir \
    --ncpu 16 2>&1 | tee ~/bipspi_run/esm2_splits/run_pathB_$(date +%Y%m%d_%H%M).log
```
Resources: **16 cpu / 192G / 12h** (codify OOM'd at 64G → needs 192G; codify does NOT cache so over-request time). No `module load` needed (ESM2 mode skips psiblast). Stages: computeFeatures (no-op) → codify(seq) → train(seq) → codify(seq_2, feedback) → train(seq_2) → eval. The **pair-AUROC table** in the log is the deliverable (vs 07b RF ~0.55). Monitor via the tee'd log + `top`/`free` over `ssh haskell` (`squeue`/`sacct` only on login node `rust`).
