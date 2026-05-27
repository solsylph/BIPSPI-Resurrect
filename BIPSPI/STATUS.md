# BIPSPI Python 3 Port — Status & Cluster Handoff

**Upstream:** [rsanchezgarc/BIPSPI](https://github.com/rsanchezgarc/BIPSPI) (Python 2.7)
**Fork purpose:** sequence-mode resurrection on the ESMFold-multimer mmseqs2-clustered splits, as a non-RF baseline after the RF-on-ESM2 architectural ceiling was confirmed (pair-AUROC ~0.55 on test).
**Last update:** 2026-05-27.

## What's been done

### 1. Mechanical Python 2 -> 3 port (complete)

21 files edited under `git diff HEAD` (+83/-63 lines). Translations:

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

### 2. Critic/implementor pair-programming round (complete)

After the mechanical port, a critic review found two correctness regressions introduced by the port:

- **P0-1 (silent wrong residue mappings):** `PairwiseAligner` in local mode returns the aligned strings TRUNCATED to the local match region; `pairwise2.align.localds` returned the full input sequences with `-` padding outside the match. Downstream code in all four alignment sites walks the strings with `idx += 1` on non-gap and indexes back into the original polypeptide — truncation shifted those indices by the leading-context offset.
- **P0-2 (crash on empty polypeptides):** `aligner.align("", X)` raises `ValueError`; `pairwise2.align.localds("", X)` returned `[]`. `homoOligomerFixer._alignSeqs` had an upstream `if len(alignments)==0: continue` guard that no longer fires.

Both fixed by an implementor pass:
- Each of the 4 alignment files (`computeFeatures/common/boundUnboundMapper.py`, `computeFeatures/common/homoOligomerFixer.py`, `pythonTools/alignSequences.py`, `computeFeatures/seqStep/seqToolManagers/seqExtraction/seqAligner.py`) gained a module-level `_padded_local_alignment(seq1, seq2, alignment)` helper that reconstructs full-length pairwise2-shaped padded strings from `alignment.aligned` block spans.
- `homoOligomerFixer._alignSeqs` gained `if not seq0 or not seq1: return []` at the top.

Byte-for-byte sanity check against legacy pairwise2 on `s1="PPPACDEFGHIKLMNQQQ", s2="ACDEFKLMN"` matched: `'PPPACDEFGHIKLMNQQQ'` / `'---ACDEF---KLMN---'`.

### 3. Cluster transfer scaffold (complete)

New files under this fork:
- `bipspi_py3_environ.yml` — Python 3.11 conda env spec (biopython >=1.81, pandas >=2.0, xgboost 1.x, gemmi, etc.).
- `tools/prepare_bipspi_inputs.py` — adapter: reads the ESMFold pipeline's `data/splits/{train,val,test}.json` + `data/structures/{PDB}_assembly1.cif`, writes BIPSPI's expected `{PREFIX}_l_b.pdb` / `{PREFIX}_r_b.pdb` per complex (plus `_u` symlinks), a single-fold `folds.json` for `--N_KFOLD <foldsFile>`, and a manifest.
- `tools/check_cluster_deps.sh` — bash probe to run on the haskell cluster login node to verify psiblast, BLAST DB, al2co, clustalw, cd-hit, SPIDER2, conda are present.
- `STATUS.md` — this file.

## What's pending

### 4. Cluster dependency probe (Gate 2)

**Not done — needs to be run on haskell.** Run:
```bash
ssh biostruct01@10.205.10.23   # via VPN-in-WSL
cd ~/BIPSPI-py3                # after you've pushed/scp'd the fork
bash tools/check_cluster_deps.sh > /tmp/bipspi_deps.txt 2>&1
cat /tmp/bipspi_deps.txt
```
Paste the output back. Sequence-mode BIPSPI requires ALL of:
- `psiblast` (NCBI BLAST+) + a compiled uniref90 BLAST DB (~60-80 GB)
- `al2co`, `clustalw`, `cd-hit` (al2co's dependencies)
- `SPIDER2` (Python 2.7 — may need its own port, or run inside a Py2 conda env subprocess)
- `conda` / `mamba` to install the Py3 env

If any are missing, that becomes a build-from-source or shim subtask before training can start.

### 5. Cluster Py3 env install

After deps pass:
```bash
mamba env create -f bipspi_py3_environ.yml      # or conda
conda activate bipspi-py3
# Smoke-import:
python -c "from Bio.Align import PairwiseAligner, substitution_matrices; print('biopython ok')"
python -c "import xgboost as xgb; print('xgboost', xgb.__version__)"
python -c "import gemmi; print('gemmi', gemmi.__version__)"
```

### 6. Smoke test on 5-10 complexes (Gate 3)

```bash
# Run the adapter against a small slice of the existing splits:
python tools/prepare_bipspi_inputs.py \
    --splits-dir ~/ESMFold-multimer/data/splits \
    --structures-dir ~/ESMFold-multimer/data/structures \
    --output-dir ~/bipspi_run/smoke \
    --evaluate test \
    --limit 5
# Then launch a tiny BIPSPI training run:
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/smoke/pdbs \
    --wdir ~/bipspi_run/smoke/wdir \
    --N_KFOLD ~/bipspi_run/smoke/folds.json \
    --ncpu 4
```
Inspect logs for: psiblast launching against uniref90, al2co running per chain, codification producing non-empty joblib pickles, xgboost training producing a model.pkl, results emitted as `prefix.tab.res.gz`.

If smoke test passes -> drop `--limit` and the full-data run is just compute time.

### 7. Full re-baseline (Gate 4)

```bash
# tmux on the login node, srun --pty inside (sbatch is blocked for regular users on haskell):
tmux new -s bipspi
srun --partition=cpu --cpus-per-task=16 --mem=64G --time=2-00:00:00 --pty bash -l
conda activate bipspi-py3
cd ~/BIPSPI-py3
python tools/prepare_bipspi_inputs.py \
    --splits-dir ~/ESMFold-multimer/data/splits \
    --structures-dir ~/ESMFold-multimer/data/structures \
    --output-dir ~/bipspi_run/full \
    --evaluate test
python generateBIPSPIModel.py \
    --modelType seq \
    --pdbsIndir ~/bipspi_run/full/pdbs \
    --wdir ~/bipspi_run/full/wdir \
    --N_KFOLD ~/bipspi_run/full/folds.json \
    --ncpu 16
# Detach: Ctrl-b d.  Reattach next day: tmux attach -t bipspi.
```

Results land under `~/bipspi_run/full/wdir/results/seq_2/` — one `.tab.res.gz` per evaluated (uppercase) prefix. Aggregate with `monitorScreenlog.py` (after fixing its pre-existing indentation bug) or by parsing the per-complex output directly.

## Reusable inputs from the ESMFold-multimer pipeline

| Pipeline artifact | Used by BIPSPI? | Through what mechanism |
|---|---|---|
| `data/splits/{train,val,test}.json` | Yes | Read directly by `tools/prepare_bipspi_inputs.py` to build the prefix list and folds JSON. |
| `data/structures/{PDB}_assembly1.cif` | Yes (with conversion) | The adapter extracts chain_a and chain_b from each mmCIF as separate PDBs using gemmi. |
| `data/labels/*_assembly1.json` | No | BIPSPI computes its own contact maps from the bound PDB during feature generation. Our JSONs aren't an input format BIPSPI knows. |
| `data/cached/esm2.zarr` | No | BIPSPI uses PSI-BLAST PSSMs, not ESM2 embeddings. |
| `data/raw/{candidates,pdb_metadata}.json` | Indirect | Useful for sequence lookup when debugging the adapter; not a direct input. |

Net saving: no re-download from RCSB.  All mmCIFs already on the cluster are reused via the adapter.

## Open decisions (deferred from Gate 1)

When you're ready to run, pick one for each:

1. **Train+val as training set, or pure train?** Default is **train+val combined as training, test held out** (matches the BIPSPI fold JSON the adapter generates). Alternative: train on train alone, evaluate on val for hyperparam tuning, then retrain train+val and evaluate on test.
2. **Metric harmonisation with 07b's RF baseline:** BIPSPI emits `auc_pair`, `prec_50/100/500`, `auc_l/r`, `mcc_l/r`. 07b emits `pair-AUROC`, `pair-AUPRC`, `prec@L/L2/L5`, `site-AUROC`, `site-F1`, `site-MCC`. Cleanest comparable: run 07b's evaluator over BIPSPI's per-complex `prefix.tab.res.gz` outputs.  Defer until after smoke test passes.

## How to transfer to the cluster

**Option A — GitHub fork:**
```bash
# locally (this Windows machine, WSL or PowerShell)
cd /e/BIPSPI-Resurrect/BIPSPI
git remote remove origin
git remote add origin git@github.com:<your-user>/BIPSPI-py3.git
git add -A && git commit -m "Python 3 port + critic fixes + cluster scaffold"
git push -u origin main
# on the cluster:
ssh biostruct01@10.205.10.23
git clone git@github.com:<your-user>/BIPSPI-py3.git ~/BIPSPI-py3
```

**Option B — direct rsync (faster for one-off, no GitHub round-trip):**
```bash
# from WSL
rsync -avz --exclude='.git/objects/pack' --exclude='__pycache__' \
    /e/BIPSPI-Resurrect/BIPSPI/ biostruct01@10.205.10.23:~/BIPSPI-py3/
```

Either works. GitHub is better if you'll iterate; rsync is one-shot. The cluster runs Linux so the relative symlinks the adapter creates for `_u.pdb -> _b.pdb` will work natively.

## Files at a glance

```
E:\BIPSPI-Resurrect\BIPSPI\
├── STATUS.md                         <- this file
├── bipspi_py3_environ.yml            <- Python 3 conda env
├── bipspi_plus_environ.yml           <- ORIGINAL Py2.7 env (kept for reference)
├── tools/
│   ├── prepare_bipspi_inputs.py      <- splits + cif -> BIPSPI pdbsIndir
│   └── check_cluster_deps.sh         <- haskell dep probe
├── Config.py                         <- PORTED
├── generateBIPSPIModel.py            <- (unchanged top-level entry)
├── predictComplexes.py               <- PORTED
├── monitorScreenlog.py               <- PORTED (with one pre-existing parse error left intact)
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
```

`git diff HEAD` shows every line of the port relative to upstream.
