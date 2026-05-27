#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 13:11:17 2020

@author: ruben
"""
import sys
from Bio.Align import PairwiseAligner, substitution_matrices

_BLOSUM62 = substitution_matrices.load("BLOSUM62")

def _padded_local_alignment(seq1, seq2, alignment):
  """Reconstruct full-length pairwise2-style padded alignment strings
  from a Bio.Align.PairwiseAligner local alignment. Returns (s1_padded, s2_padded);
  both have length len(seq1) + len(seq2) - matched_residues with '-' padding
  outside the locally aligned blocks, matching the legacy pairwise2.align.localds
  output shape."""
  t_blocks, q_blocks = alignment.aligned
  out1, out2 = [], []
  prev_t = 0
  prev_q = 0
  for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
    out1.append(seq1[prev_t:ts])
    out2.append("-" * (ts - prev_t))
    out1.append("-" * (qs - prev_q))
    out2.append(seq2[prev_q:qs])
    out1.append(seq1[ts:te])
    out2.append(seq2[qs:qe])
    prev_t, prev_q = te, qe
  out1.append(seq1[prev_t:])
  out2.append("-" * (len(seq1) - prev_t))
  out1.append("-" * (len(seq2) - prev_q))
  out2.append(seq2[prev_q:])
  return "".join(out1), "".join(out2)

def getMatchingSeqsIndices( targetSeq, referenceSeq, maxMismatchFracAllowed=0.1):
  '''
    targetSeq= "aa0aa1aa2..."  #Sequences whose indices that match referenceSeq we want to obtain
    referenceSeq= "aa0aa1aa2..."

    return: matchingIndices  target  --> ref
  '''
  targetSeq= "".join( targetSeq).replace("-","X")
  referenceSeq= "".join( referenceSeq).replace("-","X")
  matchingIndices= {}
  aligner = PairwiseAligner()
  aligner.mode = "local"
  aligner.substitution_matrix = _BLOSUM62
  aligner.open_gap_score = -11
  aligner.extend_gap_score = -0.5
  alignments = aligner.align(targetSeq, referenceSeq)
  ali_score = alignments[0].score
  target_ali, ref_ali = _padded_local_alignment(targetSeq, referenceSeq, alignments[0])
  target_ix=-1
  ref_ix=-1
  nMatches=0
  nMismatches=0
  for aa_target, aa_ref in zip(target_ali, ref_ali):
    if aa_target !="-":
      target_ix+=1
    if aa_ref !="-":
      ref_ix+=1
    if aa_target !="-" and aa_ref !="-":
      matchingIndices[target_ix]= ref_ix
      if aa_target==aa_ref:
        nMatches+=1
      else:
        nMismatches+=1

  seqIdentity= nMatches/float(nMatches+nMismatches)
  nXs= sum([1 for elem in targetSeq if elem=="X"])
  if nXs> len(targetSeq)*0.5: #skip seq if contains more than 50% of non aminoacids
    raise ValueError("To many X's in sequence")
  assert abs(len(matchingIndices) - len(referenceSeq)) < maxMismatchFracAllowed*len(referenceSeq
                                    ), "Error, excesive mismatch between \n%s\n+++++++++++++++++\n%s\n"%(targetSeq, referenceSeq)
  return matchingIndices, ali_score, seqIdentity

if __name__=="__main__":
  targetSeq, referenceSeq = sys.argv[1:]
  matchingIdxs= getMatchingSeqsIndices(targetSeq, referenceSeq)
  print(matchingIdxs)