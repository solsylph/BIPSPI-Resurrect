from __future__ import absolute_import
import os
import hashlib
import numpy as np
import pandas as pd
from ...seqToolManager import SeqToolManager
from utils import myMakeDir

class Esm2Manager(SeqToolManager):
  """
  Loads pre-computed ESM2 per-token embeddings from a Zarr store and writes
  per-chain .tab.gz feature files in BIPSPI's standard single-chain format.

  Zarr layout (written by 06_cache_esm2.py):
    store[sha256(sequence)]['per_tok_embedding']  -> float32 [L, d]
    store[sha256(sequence)]['mean_embedding']     -> float32 [d]

  Output: <computedFeatsRootDir>/esm2/<extendedPrefix>.esm2.tab.gz
    columns: chainId  resId  resName  esm2_0  esm2_1  ...  esm2_{d-1}
  """

  def __init__(self, computedFeatsRootDir, zarrPath, winSize=None):
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize)
    self.zarrPath = zarrPath
    self._store = None
    self.esm2OutPath = myMakeDir(self.computedFeatsRootDir, "esm2")

  def _getStore(self):
    if self._store is None:
      import zarr  # imported lazily so non-ESM2 runs don't require zarr
      self._store = zarr.open(self.zarrPath, mode='r')
    return self._store

  @staticmethod
  def _seqHash(seq):
    return hashlib.sha256(seq.encode()).hexdigest()

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
    seqStr = seqStr.strip()

    seq_hash = self._seqHash(seqStr)

    # Load per-token embedding from zarr store
    store = self._getStore()
    if seq_hash not in store:
      raise ValueError(
        "ESM2 embedding not found in zarr for prefix %s "
        "(hash %s...); sequence may differ from the one used during caching"
        % (prefixExtended, seq_hash[:12]))

    emb = store[seq_hash]['per_tok_embedding'][:]  # [L, d]

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
