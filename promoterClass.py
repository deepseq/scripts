########################################################################
#                                                                      #
#  This script is to classify promoters by CpG density and GC content  #
#                                                                      #
#  CpG ratio = (No. of CpGs * No. of bases)/(No. of Cs * No. of Gs)    #
#                                                                      #
########################################################################

#!/usr/bin/python3

################
#### ReadMe ####
################

'''
USAGE:\tpromoterClass.py -g <reference file> -p <promoter> -d <CpG sysmmetric file> -m <high cutoff> -n <low cutoff> -c <GC cutoff> -w <window size> -t <thread> -o <outputfile prefix>

Arguments:
\t-g, --gfile           genome reference fasta file
\t-p, --pfile           promoter regions file (must be sorted by chr ID and start position, and + separated with -)
\t-d, --dfile           CpG sysmmetric file (created by methpipe, must be sorted by chr ID and start position)
\t-m, --covcut          CpG ratio high cutoff (determine High-CpG promoters, default: 0.75)
\t-n, --percut          CpG ratio low cutoff (determine Low-CpG promoters, default: 0.48)
\t-c, --numcut          GC content cutoff (determine High-CpG promoters, default: 0.55)
\t-w, --window          window size of scanning (default: 500bp)
\t-t, --thread          number of threads to use (default: 12)
\t-o, --outputprefix    output file prefix

Flags:
\t-h, --help            help information
\t-v, --version         version information

Developer:
\tFei Sang
\tDeep Seq, Queen's Medical School
\tSchool of Life Sciences
\tUniversity of Nottingham
\tEmail: fei.sang@nottingham.ac.uk
''' 

#######################
#### main function ####
#######################

# system module
import sys, getopt
from datetime import datetime

# process module
import math
import numpy as np
import pandas as pd
import numexpr as ne
from Bio import SeqIO
from Bio.SeqUtils import GC

# multiple thread module
from tqdm import tqdm
import multiprocessing as mp
from functools import partial

def main(argvs):

  #------------#
  # parameters #
  #------------#
  genomeFile = ""
  promoterFile = ""
  cpgFile = ""
  highCut = 0.75
  lowCut = 0.48
  gcCut = 0.55
  winSize = 500
  runThread = 12
  outputPrefix = ""

  try:
    opts, args = getopt.getopt(argvs,"vhg:p:d:m:n:c:w:t:o:",["version","help","gfile=","pfile=","dfile=","highcut=","lowcut=","gccut=","window=","thread=","outputprefix="])
    if len(opts) == 0:
      print("Please use the correct arguments, for usage type -h")
      sys.exit(2)
  except(getopt.GetoptError):
    sys.exit(2)

  for opt, arg in opts:
    if opt in ("-v", "--version"):
      verbose = True
      versions()
      sys.exit()
    elif opt in ("-h", "--help"):
      usageInfo()
      sys.exit()
    elif opt in ("-g", "--gfile"):
      genomeFile = arg
    elif opt in ("-p", "--pfile"):
      promoterFile = arg
    elif opt in ("-d", "--dfile"):
      cpgFile = arg
    elif opt in ("-m", "--highcut"):
      highCut = float(arg)
    elif opt in ("-n", "--lowcut"):
      lowCut = float(arg)
    elif opt in ("-c", "--gccut"):
      gcCut = float(arg)
    elif opt in ("-w", "--window"):
      winSize = int(arg)
    elif opt in ("-t", "--thread"):
      runThread = int(arg)
    elif opt in ("-o", "--outputprefix"):
      outputPrefix = arg
    else:
      assert False, "unhandled option"

  if genomeFile == "": sys.exit("Please use the correct arguments, option -g missing")
  if promoterFile == "": sys.exit("Please use the correct arguments, option -p missing")
  if cpgFile == "": sys.exit("Please use the correct arguments, option -d missing")
  if outputPrefix == "": sys.exit("Please use the correct arguments, option -o missing")

  #------------#
  # processing #
  #------------#
  welcomeWords()
  runTime("Preparing files, please waiting ....")

  promoterMat = pd.read_csv(promoterFile, sep='\t', header=None, names=["chr","start","end","strand","geneid","genename"])
  promoterMat = promoterMat.infer_objects()

  cpgMat = pd.read_csv(cpgFile, sep='\t', header=None, names=["chr","pos","strand","tag","methfreq","basecovg"])
  cpgMat = cpgMat.infer_objects()

  runTime("Processing prometer regions by window size, please waiting ....")
  promoterMat_Win = runParallel_PromoterWin(promoterMat, runPromoterWin, winSize, runThread)

  print("Results Looking Like:")
  print(promoterMat_Win.head())

  runTime("Counting CpG sites for each region, please waiting ....")
  promoterMat_WinCpG = runParallel_Count(promoterMat_Win, runCount, cpgMat, runThread)

  print("Results Looking Like:")
  print(promoterMat_WinCpG.head())

  runTime("Calculating CpG ratio and GC content for each region, please waiting ....")
  promoterMat_WinCpG_Class = runParallel_Class(promoterMat_WinCpG, runClass, genomeFile, highCut, lowCut, gcCut, runThread)

  print("Results Looking Like:")
  print(promoterMat_WinCpG_Class.head())

  runTime("Summarizing final results, please waiting ....")
  promoterFinal = []
  regionClass = ""

  promoterRegion = promoterMat_WinCpG_Class['index'].unique().tolist()
  for idx in promoterRegion:
    regionClassMat = promoterMat_WinCpG_Class[promoterMat_WinCpG_Class['index'] == idx]
        
    if "LCP" in regionClassMat['class'].values : regionClass = "LCP"
    if "ICP" in regionClassMat['class'].values : regionClass = "ICP"
    if "HCP" in regionClassMat['class'].values : regionClass = "HCP"

    posInfo = idx.split("_")
    
    promoterFinal.append([posInfo[0],
                          posInfo[1],
                          posInfo[2],
                          regionClassMat['strand'][0],
                          regionClassMat['geneid'][0],
                          regionClassMat['genename'][0],
                          regionClass])

  promoterFinal = pd.DataFrame(promoterFinal)
  promoterFinal = promoterFinal.rename(index=str, columns={0 : "chr",
                                                           1 : "start",
                                                           2 : "end",
                                                           3 : "strand",
                                                           4 : "geneid",
                                                           5 : "genename",
                                                           6 : "class"})
  promoterFinal[['start']] = promoterFinal[['start']].astype(int)
  promoterFinal[['end']] = promoterFinal[['end']].astype(int)

  outputFile = outputPrefix + ".promoterclass.txt"
  promoterFinal.to_csv(outputFile, sep = '\t', encoding='utf-8', index=False)

  outputFile = outputPrefix + ".promoterclass.details.txt"
  promoterMat_WinCpG_Class.to_csv(outputFile, sep = '\t', encoding='utf-8', index=False)

#####################
#### subfunction ####
#####################

def versions():
  verStr = "Program:\tpromoterClass.py\nVersion:\t1.0"
  print(verStr)

def usageInfo():
  versions()
  print(__doc__)

def welcomeWords():
  welWords = "Welcome to use promoterClass.py"
  print("*"*(len(welWords)+20))
  print("*"," "*(7),welWords," "*(7),"*")
  print("*"," "*(len(welWords)+16),"*")
  print("*"," "*(12),"Fei Sang,DeepSeq,QMC"," "*(13),"*")
  print("*"," "*(6),"Life Sciences,Univ. of Nottingham"," "*(6),"*")
  print("*"," "*(11),"Nottingham, NG7 2UH, UK"," "*(11),"*")
  print("*"," "*(10),"fei.sang@nottingham.ac.uk"," "*(10),"*")
  print("*"*(len(welWords)+20))

def runTime(strTmp):
  runningTime = str(datetime.now())
  runningTime = runningTime.split(".")[0]
  print(runningTime,">>>>",strTmp)

#----------------------------------------------------------------------------#
# for loop --> range(0,3), loop from 0 to 2, 3 times                         #
# sequence[0:2], index 2 base not included                                   #
# pd dataframe append too slow, use list instead, and convert into dataframe #
#----------------------------------------------------------------------------#

def runPromoterWin(promoterMat, winSize):
  promoterWin = []

  pbar = tqdm(range(len(promoterMat.index)), ncols = 100)
  for i in pbar:
    for j in range(promoterMat.iloc[i,1],(promoterMat.iloc[i,2]-winSize+1+1)):
      promoterWin.append([promoterMat.iloc[i,0]+"_"+str(promoterMat.iloc[i,1])+"_"+str(promoterMat.iloc[i,2]),
                          promoterMat.iloc[i,0],
                          j,
                          j+winSize-1,
                          promoterMat.iloc[i,3],
                          promoterMat.iloc[i,4],
                          promoterMat.iloc[i,5]])
  pbar.close()

  promoterWin = pd.DataFrame(promoterWin)
  promoterWin = promoterWin.rename(index=str, columns={0 : "index",
                                                       1 : "chr",
                                                       2 : "start",
                                                       3 : "end",
                                                       4 : "strand",
                                                       5 : "geneid",
                                                       6 : "genename"})
  promoterWin[['start']] = promoterWin[['start']].astype(int)
  promoterWin[['end']] = promoterWin[['end']].astype(int)
  return promoterWin

def runCount(promoterMat_Win, cpgMat):
  promoterCpG = []
  tmpchr = ""

  pbar = tqdm(range(len(promoterMat_Win.index)), ncols = 100)
  for i in pbar:
    if tmpchr != promoterMat_Win.iloc[i,1] :
      cpgMat_tmp = cpgMat.loc[cpgMat['chr'] == promoterMat_Win.iloc[i,1]]
      tmpchr = promoterMat_Win.iloc[i,1]
    
    cpgMat_region_values = cpgMat_tmp['pos'].values
    v1 = promoterMat_Win.iloc[i,2]
    v2 = promoterMat_Win.iloc[i,3]
    cpgMat_region = cpgMat_tmp[ne.evaluate("(cpgMat_region_values >= v1) & (cpgMat_region_values <= v2)")].values

    #cpgMat_region = cpgMat_tmp.loc[(cpgMat_tmp['pos'] >= promoterMat_Win.iloc[i,2]) &
    #                               (cpgMat_tmp['pos'] <= promoterMat_Win.iloc[i,3])].values

    cpgCount = len(cpgMat_region)

    promoterCpG.append([promoterMat_Win.iloc[i,0],
                        promoterMat_Win.iloc[i,1],
                        promoterMat_Win.iloc[i,2],
                        promoterMat_Win.iloc[i,3],
                        promoterMat_Win.iloc[i,4],
                        promoterMat_Win.iloc[i,5],
                        promoterMat_Win.iloc[i,6],
                        cpgCount])
  pbar.close()

  promoterCpG = pd.DataFrame(promoterCpG)
  promoterCpG = promoterCpG.rename(index=str, columns={0 : "index",
                                                       1 : "chr",
                                                       2 : "start",
                                                       3 : "end",
                                                       4 : "strand",
                                                       5 : "geneid",
                                                       6 : "genename",
                                                       7 : "cpgcount"})
  promoterCpG[['start']] = promoterCpG[['start']].astype(int)
  promoterCpG[['end']] = promoterCpG[['end']].astype(int)
  promoterCpG[['cpgcount']] = promoterCpG[['cpgcount']].astype(int)
  return promoterCpG

def runClass(promoterMat_WinCpG, genomeFile, highCut, lowCut, gcCut):
  genomeDict = SeqIO.to_dict(SeqIO.parse(genomeFile, "fasta"))
  promoterClass = []

  pbar = tqdm(range(len(promoterMat_WinCpG.index)), ncols = 100)
  for i in pbar:
    regionSeq = genomeDict[promoterMat_WinCpG.iloc[i,1]].seq[promoterMat_WinCpG.iloc[i,2]:(promoterMat_WinCpG.iloc[i,3]+1)]
    if (regionSeq.count("C") == 0) or (regionSeq.count("G") == 0):
      cpgRatio = -0.1
    else:
      cpgRatio = (promoterMat_WinCpG.iloc[i,7] * len(regionSeq))/(regionSeq.count("C") * regionSeq.count("G"))
    gcContent = GC(regionSeq)/100

    if (cpgRatio >= highCut) & (gcContent >= gcCut): regionClass = "HCP"
    elif cpgRatio < lowCut : regionClass = "LCP"
    else : regionClass = "ICP"
  
    promoterClass.append([promoterMat_WinCpG.iloc[i,0],
                          promoterMat_WinCpG.iloc[i,1],
                          promoterMat_WinCpG.iloc[i,2],
                          promoterMat_WinCpG.iloc[i,3],
                          promoterMat_WinCpG.iloc[i,4],
                          promoterMat_WinCpG.iloc[i,5],
                          promoterMat_WinCpG.iloc[i,6],
                          promoterMat_WinCpG.iloc[i,7],
                          cpgRatio,
                          gcContent,
                          regionClass])
  pbar.close()
  
  promoterClass = pd.DataFrame(promoterClass)
  promoterClass = promoterClass.rename(index=str, columns={0  : "index",
                                                           1  : "chr",
                                                           2  : "start",
                                                           3  : "end",
                                                           4  : "strand",
                                                           5  : "geneid",
                                                           6  : "genename",
                                                           7  : "cpgcount",
                                                           8  : "cpgratio",
                                                           9  : "gc",
                                                           10 : "class"})
  promoterClass[['start']] = promoterClass[['start']].astype(int)
  promoterClass[['end']] = promoterClass[['end']].astype(int)
  promoterClass[['cpgcount']] = promoterClass[['cpgcount']].astype(int)
  promoterClass[['cpgratio']] = promoterClass[['cpgratio']].astype(float)
  promoterClass[['gc']] = promoterClass[['gc']].astype(float)
  return promoterClass

def runParallel_PromoterWin(df, func, w, t):
  df_split = np.array_split(df, t)
  pool = mp.Pool(t)
  df_out = pool.map(partial(func, winSize=w), df_split)
  pool.close()
  pool.join()
  return pd.concat(df_out)

def runParallel_Count(df, func, d, t):
  df_split = np.array_split(df, t)
  pool = mp.Pool(t)
  df_out = pool.map(partial(func, cpgMat=d), df_split)
  pool.close()
  pool.join()
  return pd.concat(df_out)

def runParallel_Class(df, func, g, h, l, c, t):
  df_split = np.array_split(df, t)
  pool = mp.Pool(t)
  df_out = pool.map(partial(func, genomeFile=g, highCut=h, lowCut=l, gcCut=c), df_split)
  pool.close()
  pool.join()
  return pd.concat(df_out)

#####################
#### program run ####
#####################

if __name__ == "__main__":
  main(sys.argv[1:])

