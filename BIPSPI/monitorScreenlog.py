import sys, os
import pandas as pd

pd.set_option('display.precision', 4)

EVAL_JUST_LAST_STEP=False

SKIP_LOWERCASE=True

def main(fname):
  pdbsScoresDict={}
  nextLineIsGood=False
  colNames= None
  with open(fname) as f:
    for line in f:
#      print(line)
      if "pdb  auc_pair  " in line:
        colNames= line.split()
        nextLineIsGood=True
        continue
      if nextLineIsGood:
        dataLine= line.split()
        pdbId= dataLine[0].split(".")[0]
        if EVAL_JUST_LAST_STEP and "@" in pdbId:
          nextLineIsGood=False
          continue
        if SKIP_LOWERCASE and pdbId.islower(): continue
          nextLineIsGood=False
          continue
        dataLine= [dataLine[0]]+ [float(elem) for elem in dataLine[1:]]
        pdbsScoresDict[dataLine[0]]= dataLine
        nextLineIsGood=False
  resDf= pd.DataFrame.from_records(list(pdbsScoresDict.values()))
  resDf.columns= colNames
  resDf.sort_values("pdb", inplace=True)

  means= resDf.mean(axis=0)
  resDf= pd.concat([resDf, resDf.iloc[[-1]]], ignore_index=True)
  resDf.iloc[-1, 0]=  "mean"
  resDf.iloc[-1, 1:]=  means
  print(resDf.to_string(index=False))
  
if __name__=="__main__":
  if len(sys.argv)==2:
    fname=sys.argv[1]
  else:
    fname= "screenlog.0"
  main(fname)
