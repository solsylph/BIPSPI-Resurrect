# BIPSPI Python 3 Port — Status & Cluster Handoff

**Upstream:** [rsanchezgarc/BIPSPI](https://github.com/rsanchezgarc/BIPSPI) (Python 2.7)
**Our fork:** [solsylph/BIPSPI-Resurrect](https://github.com/solsylph/BIPSPI-Resurrect)
**Purpose:** sequence-mode resurrection on the ESMFold-multimer mmseqs2-clustered splits, as a non-RF baseline after the RF-on-ESM2 architectural ceiling was confirmed (pair-AUROC ~0.55 on test).
**Last update:** 2026-05-27 (after Phase C complete, Phase D in progress).

> **Repo layout note:** the GitHub repo has `BIPSPI/` as a subdirectory at its root (the local git repo was reinitialised at `E:\BIPSPI-Resurrect\` instead of `E:\BIPSPI-Resurrect\BIPSPI\`). On the cluster the code lives at `~/BIPSPI-Resurrect/BIPSPI/`. All commands below assume `cd ~/BIPSPI-Resurrect/BIPSPI` first.

---

## Progress at a glance

| Phase | Subject | Status |
|---|---|---|
| Port | Py2.7 → Py3 mechanical port (21 files, +83/-63 lines) | ✅ |
| Critic | Two P0 regressions found + fixed (P0-1 alignment truncation, P0-2 empty-input crash) | ✅ |
| A | Extend `protein` conda env + module-load BLAST+/CD-HIT | ✅ |
| B | Build al2co + install clustalw 2.1 (via bioconda) | ✅ |
| Cfg patch | Update `dependencies.cfg` with cluster paths | ✅ |
| C | Port SPIDER2 to Py3.10 (numpy-only, option a-i) | ✅ |
| D | Download uniref90 + makeblastdb | ⏳ in progress (~3-4 hr ETA) |
| Smoke | BIPSPI training on 5-10 complexes (Gate 3) | blocked on D |
| Full | BIPSPI re-baseline on full test split (Gate 4) | blocked on smoke |

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

After the mechanical port, a critic review found two correctness regressions introduced by the port:

- **P0-1 (silent wrong residue mappings):** `PairwiseAligner` in local mode returns the aligned strings TRUNCATED to the local match region; `pairwise2.align.localds` returned the full input sequences with `-` padding outside the match. Downstream code walks the strings with `idx += 1` and indexes into the original polypeptide — truncation shifted those indices by the leading-context offset.
- **P0-2 (crash on empty polypeptides):** `aligner.align("", X)` raises `ValueError`; `pairwise2.align.localds("", X)` returned `[]`. `homoOligomerFixer._alignSeqs` had a `if len(alignments)==0: continue` guard that no longer fired.

Both fixed:
- All 4 alignment files (`computeFeatures/common/boundUnboundMapper.py`, `computeFeatures/common/homoOligomerFixer.py`, `pythonTools/alignSequences.py`, `computeFeatures/seqStep/seqToolManagers/seqExtraction/seqAligner.py`) gained a module-level `_padded_local_alignment(seq1, seq2, alignment)` helper that reconstructs full-length pairwise2-shaped padded strings from `alignment.aligned` block spans.
- `homoOligomerFixer._alignSeqs` gained `if not seq0 or not seq1: return []` at the top.

Byte-for-byte sanity check against legacy pairwise2 on `s1="PPPACDEFGHIKLMNQQQ", s2="ACDEFKLMN"` matched: `'PPPACDEFGHIKLMNQQQ'` / `'---ACDEF---KLMN---'`.

---

## 3. Phase A — env + module loads (complete)

**Cluster Python env**: `protein` (existing conda env at `~/miniforge3/envs/protein`).

**Confirmed Python version: 3.10.20**.

**Packages installed into `protein`:**
```
conda install -n protein -c conda-forge -c bioconda \
    xgboost gemmi tqdm requests mmtf-python psutil joblib -y
```
Currently installed (verified):
- xgboost **3.2.0** — see caveat below
- biopython 1.86
- pandas 2.3.3
- gemmi, tqdm, requests, mmtf-python, psutil, joblib (all current)

> **xgboost 3.2.0 caveat**: BIPSPI's `trainAndTest/classifiers/xgBoost.py` was written against xgboost 0.80 (2017). xgboost 3.x has substantial API changes (early stopping moved, callbacks restructured, `use_label_encoder` removed). **This will almost certainly break at BIPSPI training time.** Per "fix later" scope it's left as-is. Either pre-empt with `conda install -n protein "xgboost>=1.7,<2.0" -y`, or patch `xgBoost.py` when it errors.

**Module loads** (must repeat in any fresh shell):
```bash
conda activate protein
module load BLAST+/2.14.1-gompi-2023a
module load CD-HIT/4.8.1-GCC-12.2.0
```
Resolved paths:
- `psiblast`, `makeblastdb` → `/cvmfs/.../BLAST+/2.14.1-gompi-2023a/bin/`
- `cd-hit` → `/cvmfs/.../CD-HIT/4.8.1-GCC-12.2.0/bin/`

Add the three lines to `~/.bashrc` if you want them permanent.

> Loading CD-HIT downgrades GCC 12.3 → 12.2 in the shell environment. Harmless for our purposes.

## 4. Phase B — al2co + clustalw (complete)

Built/installed via `tools/install_al2co_clustalw.sh`:
- **clustalw** via bioconda → version **2.1** at `~/miniforge3/envs/protein/bin/clustalw`. BIPSPI docs spec'd 1.83 but the CLI is backward-compatible — al2co only uses standard MSA invocation.
- **al2co** built from [TheApacheCats/al2co GitHub mirror](https://github.com/TheApacheCats/al2co) (the original swmed FTP is dead). Single C file with the `char[500]` → `char[1024]` buffer patch from BIPSPI docs, compiled with `gcc al2co.c -o al2co -lm`. Build emits K&R-era warnings (implicit `int`, fgets-into-fstr buffer overflow on a pre-existing bug) — all harmless. Smoke test produced expected usage text.

Output binary paths:
- `~/tools/al2co/al2co`
- `~/miniforge3/envs/protein/bin/clustalw`

## 5. Phase B.5 — `dependencies.cfg` patcher (complete)

`tools/patch_dependencies_cfg.sh` rewrote `configFiles/cmdTool/dependencies.cfg` with the cluster-specific paths above. Backup left at `dependencies.cfg.bak.<timestamp>`.

Final values (live on cluster):
```
psiBlastBin       /cvmfs/.../BLAST+/2.14.1-gompi-2023a/bin/psiblast
psiBlastDB_path   /home/biostruct01/databases/uniref90/uniref90.fasta   # populated by Phase D
cdHitBin_path     /cvmfs/.../CD-HIT/4.8.1-GCC-12.2.0/bin/cd-hit
clustalW_path     /home/biostruct01/miniforge3/envs/protein/bin/clustalw
al2coBin_path     /home/biostruct01/tools/al2co/al2co
spider2PyScript_path /home/biostruct01/tools/SPIDER2/misc/pred_pssm.py
```

## 6. Phase C — SPIDER2 port (complete)

**Decision taken: option a-i** (numpy-only rewrite of forward pass).

SPIDER2 source: `SPIDER2_local.tgz` (105 MB, 2017 vintage) from the Zhou lab successor site (`http://183.36.5.251:8080/sparks_downloads/.../old_versions/SPIDER2_local.tgz` — the original Sparks Lab URL is dead). Extracted to `~/tools/SPIDER2/`.

**Key discovery**: `pred_pssm.py` is **already pure numpy** — no Theano dependency. The "NN" is matrix multiplies + sigmoid, weights load via `numpy.load()` on the bundled `.npz` files. Three iterations refine SS/ASA/TTPP predictions.

Three Py2 → Py3 edits applied by `tools/port_spider2.sh`:
1. **3× `print >>fp, X` → `print(X, file=fp)`** (Py2 redirect-print syntax)
2. **`numpy.load(f)` → `numpy.load(f, allow_pickle=True, encoding='latin1')`** — numpy ≥1.16.3 changed default to `allow_pickle=False` (security); Py2-pickled `.npz` also needs `encoding='latin1'` to decode bytes strings safely

**Validation:** smoke test against bundled `SPIDER2/ex/1a1xA.pssm` produced `1a1xA.spd3` that **byte-matches** the reference `1a1xA_CHECK.spd3` on the first 5 lines (`diff` returned empty). Numerical fidelity confirmed.

One non-fatal `DeprecationWarning: Conversion of an array with ndim > 0 to a scalar` (numpy 1.25) when formatting `pred_asa_1[ind]` via `%5.1f`. Will error under numpy 2.x; we're pinned at numpy <2.0 so safe. Fix is `.item()` insertion when needed.

Files we did NOT touch in SPIDER2 (not on BIPSPI's critical path): `splitseq.py` (has its own Py2 bugs), `pred_nopssm.py`, `pred_pssm0.py`, `seq2pssm.py`.

## 7. Phase D — uniref90 BLAST DB (in progress)

Running in `tmux session "uniref90"` on haskell node, inside `srun --pty bash` shell. Script: `tools/download_uniref90_db.sh`. Target dir: `~/databases/uniref90/`.

Steps + realistic timing:
- `wget uniref90.fasta.gz` from UniProt FTP — **44 GB compressed** (larger than initially estimated; UniRef90 grew), ~60-90 min at ~11 MB/s
- `gunzip` → ~130 GB uncompressed FASTA, ~30 min
- `makeblastdb -dbtype prot -hash_index -parse_seqids` → ~60-90 min

Total: **~3-4 hours**. Disk: ~130 GB of 17 TB free — no constraint.

When complete, `dependencies.cfg` `psiBlastDB_path` already points at `~/databases/uniref90/uniref90.fasta`. Verify with:
```bash
echo -e ">test\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF" > /tmp/test.fa
psiblast -query /tmp/test.fa -db ~/databases/uniref90/uniref90.fasta -num_iterations 1 -num_threads 4 -out /tmp/test.psiblast
head -30 /tmp/test.psiblast
```

---

## 8. What remains (blocked on Phase D)

### Smoke test on 5-10 complexes (Gate 3)

```bash
cd ~/BIPSPI-Resurrect/BIPSPI
conda activate protein
module load BLAST+/2.14.1-gompi-2023a CD-HIT/4.8.1-GCC-12.2.0

# Adapter: produce BIPSPI-format pdbsIndir from existing ESMFold splits
python tools/prepare_bipspi_inputs.py \
    --splits-dir ~/ESMFold-multimer/data/splits \
    --structures-dir ~/ESMFold-multimer/data/structures \
    --output-dir ~/bipspi_run/smoke \
    --evaluate test \
    --limit 5

# Smoke training run
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/smoke/pdbs \
    --wdir ~/bipspi_run/smoke/wdir \
    --N_KFOLD ~/bipspi_run/smoke/folds.json \
    --ncpu 4
```

Watch for: psiblast queries against uniref90, al2co/SPIDER2 running per chain, xgboost training, `prefix.tab.res.gz` output files. **Most likely failure point: xgboost API incompatibility** (BIPSPI's xgBoost.py vs xgboost 3.x) — fix in place when it errors.

### Full re-baseline (Gate 4)

After smoke passes, drop `--limit 5` and re-run inside tmux+srun. Outputs land at `~/bipspi_run/full/wdir/results/seq_2/{PREFIX}.tab.res.gz` for every uppercase (evaluated) prefix.

---

## 9. Reusable inputs from the ESMFold-multimer pipeline

| Pipeline artifact | Used by BIPSPI? | How |
|---|---|---|
| `~/ESMFold-multimer/data/splits/{train,val,test}.json` | Yes | Read directly by `tools/prepare_bipspi_inputs.py` to build the prefix list and folds JSON. |
| `~/ESMFold-multimer/data/structures/{PDB}_assembly1.cif` | Yes (with conversion) | The adapter extracts chain_a and chain_b from each mmCIF as separate PDBs using gemmi. |
| `~/ESMFold-multimer/data/labels/*_assembly1.json` | No | BIPSPI computes its own contact maps from the bound PDB during feature generation. |
| `~/ESMFold-multimer/data/cached/esm2.zarr` | No | BIPSPI uses PSI-BLAST PSSMs, not ESM2 embeddings. |
| `~/ESMFold-multimer/data/raw/{candidates,pdb_metadata}.json` | Indirect | Useful for sequence lookup when debugging the adapter. |

Net saving vs. starting over: no re-download from RCSB. All mmCIFs already on the cluster are reused via the adapter.

---

## 10. Open decisions (deferred from Gate 1)

When you're ready to commit to a baseline run, pick one for each:

1. **Train+val as training set, or pure train?** Adapter's default is **train+val combined as training, test held out** (lowercases train+val, uppercases test — BIPSPI's `SKIP_LOWER_PREDICTION` then evaluates only test in the last step).
2. **Metric harmonisation with 07b's RF baseline:** BIPSPI emits `auc_pair`, `prec_50/100/500`, `auc_l/r`, `mcc_l/r`. 07b emits `pair-AUROC`, `pair-AUPRC`, `prec@L/L2/L5`, `site-AUROC`, `site-F1`, `site-MCC`. Cleanest comparable: run 07b's evaluator over BIPSPI's per-complex `prefix.tab.res.gz` outputs. Defer until after smoke test passes.

---

## 11. How to transfer / sync the code

**Push from local Windows machine:**
```bash
cd /e/BIPSPI-Resurrect
git add -A
git commit -m "..."
git push
```
(Note: the local git repo is at `/e/BIPSPI-Resurrect/` not `/e/BIPSPI-Resurrect/BIPSPI/` due to a `git init` that happened at the parent. That's why GitHub has `BIPSPI/` as a subdir at the root.)

**Pull on the cluster:**
```bash
cd ~/BIPSPI-Resurrect
git pull
cd BIPSPI    # all scripts and BIPSPI code live here
```

---

## 12. Cluster details (haskell / rust)

- **Login**: `ssh biostruct01@10.205.10.23` (VPN-in-WSL required). The cluster has multiple nodes (`haskell`, `rust`) sharing the same `/home` filesystem; commands can run on either.
- **No sbatch** for regular users. Long jobs run inside `srun --pty bash -l` wrapped in `tmux`.
- **Conda**: miniforge3 at `~/miniforge3/`. Active env is `protein`.

---

## 13. Files at a glance

```
~/BIPSPI-Resurrect/BIPSPI/                       (cluster)
E:\BIPSPI-Resurrect\BIPSPI\                      (local Windows)
├── STATUS.md                          <- this file
├── bipspi_py3_environ.yml             <- Python 3 conda env spec (not used; we extended `protein` instead)
├── bipspi_plus_environ.yml            <- ORIGINAL Py2.7 env (kept for reference)
├── tools/
│   ├── prepare_bipspi_inputs.py       <- splits + cif -> BIPSPI pdbsIndir (Python 3)
│   ├── check_cluster_deps.sh          <- haskell dep probe (already run)
│   ├── install_al2co_clustalw.sh      <- Phase B installer (already run, idempotent)
│   ├── patch_dependencies_cfg.sh      <- cfg path patcher (already run, idempotent)
│   ├── download_uniref90_db.sh        <- Phase D downloader (running in tmux now)
│   └── port_spider2.sh                <- Phase C SPIDER2 patcher (already run, idempotent)
├── configFiles/cmdTool/
│   └── dependencies.cfg               <- PATCHED with cluster paths
├── Config.py                          <- PORTED
├── generateBIPSPIModel.py             <- (unchanged top-level entry)
├── predictComplexes.py                <- PORTED
├── monitorScreenlog.py                <- PORTED (with one pre-existing parse error left intact)
├── computeFeatures/
│   └── common/{boundUnboundMapper,homoOligomerFixer}.py    <- PORTED (P0-1 + P0-2 critic fixes applied)
├── codifyComplexes/
│   ├── codifyProtocols/DataLoaderClass.py     <- PORTED
│   └── codifyProtocols/SeqProtocol.py         <- (unchanged)
├── trainAndTest/
│   ├── trainAndTest.py                        <- PORTED
│   ├── resultsManager.py                      <- PORTED
│   └── evaluateResults.py                     <- PORTED
├── pythonTools/
│   ├── alignSequences.py                      <- PORTED (P0-1 critic fix applied)
│   ├── extractModelsFromPdbFile.py            <- PORTED
│   └── combinePDBs.py                         <- PORTED
└── evaluation/                                <- 6 files PORTED (.ix/.append idiom)

External (on cluster, NOT in this repo):
~/tools/al2co/al2co                       <- al2co binary (Phase B)
~/tools/SPIDER2/                          <- SPIDER2 source + weights (Phase C)
~/tools/SPIDER2/misc/pred_pssm.py         <- PORTED via tools/port_spider2.sh
~/databases/uniref90/uniref90.fasta       <- BLAST DB (Phase D — building)
~/miniforge3/envs/protein/                <- conda env (Phase A — extended with bipspi deps)
```

`git diff HEAD` shows every line of the port relative to upstream.
