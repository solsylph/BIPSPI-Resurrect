from __future__ import absolute_import
import os
import json
import hashlib
import numpy as np
import pandas as pd
from ...seqToolManager import SeqToolManager
from utils import myMakeDir

# Minimum fraction of a chain's residues that must map onto a candidate
# embedding for the SEQRES fallback to be accepted.
MIN_COVERAGE = 0.95

class Esm2Manager(SeqToolManager):
  """
  Loads pre-computed ESM2 per-token embeddings from a Zarr store and writes
  per-chain .tab.gz feature files in BIPSPI's standard single-chain format.

  Zarr layout (written by 06_cache_esm2.py):
    store[sha256(sequence)[:16]]['per_tok_embedding']  -> float32 [L, d]
    store[sha256(sequence)[:16]]['mean_embedding']     -> float32 [d]

  The zarr is keyed on the SHA256 (first 16 hex chars) of the SEQRES sequence
  used when the cache was built (from candidates.json). BIPSPI, however,
  extracts sequences from PDB ATOM records, which omit unresolved residues, so
  the ATOM-sequence hash often misses. When it does and a candidates.json is
  provided, we look up the SEQRES sequence for the same pdbId, fetch its
  per-token embedding, and align+slice it down to BIPSPI's resolved residues.

  Output: <computedFeatsRootDir>/esm2/<extendedPrefix>.esm2.tab.gz
    columns: chainId  resId  resName  esm2_0  esm2_1  ...  esm2_{d-1}
  """

  def __init__(self, computedFeatsRootDir, zarrPath, candidatesJsonPath=None, winSize=None):
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)
    self.zarrPath = zarrPath
    self.candidatesJsonPath = candidatesJsonPath
    self._store = None
    self._candidatesByPdb = None  # lazily built {pdbIdLower: [seqresSeq, ...]}
    self.esm2OutPath = myMakeDir(self.computedFeatsRootDir, "esm2")

  def _getStore(self):
    if self._store is None:
      import zarr  # imported lazily so non-ESM2 runs don't require zarr
      self._store = zarr.open(self.zarrPath, mode='r')
    return self._store

  @staticmethod
  def _seqHash(seq):
    # Must match the key derivation in 06_cache_esm2.py (the zarr writer),
    # which truncates the SHA256 hexdigest to 16 chars.
    return hashlib.sha256(seq.encode()).hexdigest()[:16]

  def _getCandidates(self):
    '''
      Lazily parse candidates.json into {pdbIdLower: [seqresSeq, ...]}.
      Returns {} if no candidates file is configured.
    '''
    if self._candidatesByPdb is not None:
      return self._candidatesByPdb
    self._candidatesByPdb = {}
    if not self.candidatesJsonPath or not os.path.isfile(self.candidatesJsonPath):
      return self._candidatesByPdb
    with open(self.candidatesJsonPath) as f:
      data = json.load(f)
    cands = data.get('candidates', data if isinstance(data, list) else [])
    for entry in cands:
      pdbId = str(entry.get('pdb_id', '')).lower()
      if not pdbId:
        continue
      seqs = self._candidatesByPdb.setdefault(pdbId, [])
      for chain in entry.get('chains', []):
        seq = chain.get('sequence', '')
        if seq:
          seqs.append(seq)
    return self._candidatesByPdb

  @staticmethod
  def _alignAndSlice(atomSeq, seqresSeq, perTokEmb, meanEmb):
    '''
      Align the (resolved) ATOM sequence to the (full) SEQRES sequence and slice
      the SEQRES per-token embedding down to one row per ATOM residue, in order.
      Unmapped ATOM residues (rare; ATOM not present in SEQRES) are filled with
      the chain mean embedding.

      :return (alignedEmb [len(atomSeq), d], nMapped) or (None, 0) if alignment
              is too poor to trust.
    '''
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    # Identity-style scoring: robust to any residue letter (X, U, B, ...) that a
    # substitution matrix might not contain. ATOM is (almost) a subsequence of
    # SEQRES, so we just need to recover the position mapping.
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    try:
      aln = aligner.align(seqresSeq, atomSeq)[0]  # target=SEQRES, query=ATOM
    except Exception:
      return None, 0

    d = perTokEmb.shape[1]
    alignedEmb = np.empty((len(atomSeq), d), dtype=perTokEmb.dtype)
    mapped = np.zeros(len(atomSeq), dtype=bool)

    targetBlocks, queryBlocks = aln.aligned  # equal-length aligned runs
    for (ts, te), (qs, qe) in zip(targetBlocks, queryBlocks):
      for off in range(te - ts):
        atomPos = qs + off
        seqresPos = ts + off
        alignedEmb[atomPos] = perTokEmb[seqresPos]
        mapped[atomPos] = True

    nMapped = int(mapped.sum())
    if nMapped < MIN_COVERAGE * len(atomSeq):
      return None, nMapped

    if nMapped < len(atomSeq):  # fill the few unmapped residues with the mean
      alignedEmb[~mapped] = meanEmb
    return alignedEmb, nMapped

  def _lookupEmbedding(self, prefix, seqStr, prefixExtended):
    '''
      Resolve a [len(seqStr), d] per-residue embedding for this chain, trying an
      exact ATOM-sequence hash first, then the SEQRES candidates fallback.
    '''
    store = self._getStore()

    # Fast path: ATOM sequence was cached verbatim.
    seqHash = self._seqHash(seqStr)
    if seqHash in store:
      emb = store[seqHash]['per_tok_embedding'][:]
      if emb.shape[0] == len(seqStr):
        return emb
      # hash collision-free but length disagrees: fall through to alignment

    # Fallback: align the candidate SEQRES embedding onto the ATOM residues.
    pdbId = prefix.split("@")[0].lower()
    candidateSeqs = self._getCandidates().get(pdbId, [])
    best = None  # (nMapped, alignedEmb)
    seenHashes = set()
    for cseq in candidateSeqs:
      ch = self._seqHash(cseq)
      if ch in seenHashes or ch not in store:
        continue
      seenHashes.add(ch)
      group = store[ch]
      aligned, nMapped = self._alignAndSlice(
        seqStr, cseq, group['per_tok_embedding'][:], group['mean_embedding'][:])
      if aligned is not None and (best is None or nMapped > best[0]):
        best = (nMapped, aligned)

    if best is not None:
      return best[1]

    raise ValueError(
      "ESM2 embedding not found in zarr for prefix %s "
      "(hash %s...) and no candidate SEQRES sequence for pdbId %s aligned; "
      "sequence may differ from the one used during caching"
      % (prefixExtended, seqHash[:12], pdbId))

  def getFinalPath(self):
    return self.esm2OutPath

  def getFNames(self, prefixExtended):
    return [os.path.join(self.esm2OutPath, prefixExtended + ".esm2.tab.gz")]

  def compressRawData(self, prefixExtended):
    pass  # already gzipped on write

  def computeFromSeqStructMapper(self, seqStructMap, prefixExtended):
    """
    Load ESM2 embedding for the chain identified by prefixExtended and write
    a tab-separated feature file.  Returns the output filename.
    """
    outFname = self.getFNames(prefixExtended)[0]
    if os.path.isfile(outFname):
      return outFname

    prefix, chainType, chainId = self.splitExtendedPrefix(prefixExtended)[:3]
    seqStr, _ = seqStructMap.getSeq(chainType, chainId)
    # Mirror every classic single-chain manager (PsiBlast/Spider2/Al2co/HHblits/
    # windowSeq): setCurrentSeq() registers this chain in seqStructMap.seqToRefSeq.
    # Without it, seqToStructIndex() hits `seqToRefSeq[(chainType,chainId)]` ->
    # KeyError -> returns None for EVERY residue, so resId falls back to the bogus
    # "<seqIx>?" form and the codify inner-join on resId yields 0 rows ("dataset
    # is empty"). Must be called with the seq straight from getSeq so it matches
    # the stored seqsDict entry (identity mapping).
    seqStructMap.setCurrentSeq(seqStr, chainType, chainId)
    seqStr = seqStr.strip()

    emb = self._lookupEmbedding(prefix, seqStr, prefixExtended)  # [len(seqStr), d]

    if emb.shape[0] != len(seqStr):
      raise ValueError(
        "ESM2 embedding length (%d) != sequence length (%d) for %s"
        % (emb.shape[0], len(seqStr), prefixExtended))

    d = emb.shape[1]
    col_names = ["esm2_%d" % i for i in range(d)]

    rows = []
    for seqIx in range(len(seqStr)):
      structIndex = seqStructMap.seqToStructIndex(chainType, chainId, seqIx, asString=True)
      if structIndex is None:
        structIndex = str(seqIx) + "?"
      rows.append([chainId, structIndex, seqStr[seqIx]] + emb[seqIx].tolist())

    df = pd.DataFrame(rows, columns=["chainId", "resId", "resName"] + col_names)
    df.to_csv(outFname, sep='\t', index=False, compression='gzip')
    return outFname
