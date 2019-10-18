#################################################
#                                               #
#  This script is for bin methylation levels.	#
#                                               #
#################################################

#!/usr/bin/python3

################
#### ReadMe ####
################

'''
USAGE:\tbinMethLevel.py -i <inputfile> -b <binsize> -c <covcut> -t <thread> -o <outputfile prefix>

Arguments:
\t-i, --input           input file
\t-b, --binsize         bin size
\t-c, --covcut          coverage cutoff
\t-t, --thread          number of threads to use
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

import sys, getopt
from datetime import datetime

import math
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from functools import partial

def main(argvs):

  #------------#
  # parameters #
  #------------#
  inputFile = ''
  outputPrefix = ''

  try:
    opts, args = getopt.getopt(argvs,"vhi:b:c:t:o:",["version","help","input=","binsize=","covcut=","thread=","outputprefix="])
    if len(opts) == 0:
      print("Please use the correct arguments, for usage type -h")
      sys.exit(2)
  except(getopt.GetoptError, err):
    print(str(err))
    sys.exit(2)

  for opt, arg in opts:
    if opt in ("-v", "--version"):
      verbose = True
      versions()
      sys.exit()
    elif opt in ("-h", "--help"):
      usageInfo()
      sys.exit()
    elif opt in ("-i", "--input"):
      inputFile = arg
    elif opt in ("-b", "--binsize"):
      binSize = int(arg)
    elif opt in ("-c", "--covcut"):
      covCut = int(arg)
    elif opt in ("-t", "--thread"):
      runThread = int(arg)
    elif opt in ("-o", "--outputprefix"):
      outputPrefix = arg
    else:
      assert False, "unhandled option"

  if inputFile == "": sys.exit("Please use the correct arguments, option -i missing")
  if outputPrefix == "": sys.exit("Please use the correct arguments, option -o missing")

  num_partitions = 4 #number of partitions to split dataframe
  num_cores = 4 #number of cores on your machine

  #------------#
  # processing #
  #------------#
  welcomeWords()
  runTime("Fetch CpG Methylation Level, Please Waiting ....")

  methLevelMat = pd.read_csv(inputFile, sep='\t', header=None, names=["chr","pos","strand","tag","freq","covg"])
  methLevelMat = methLevelMat.infer_objects()
  methLevelMat = methLevelMat.round({'freq':3})

  runTime("Creating CpG Methylation Level Bin Matrix, Please Waiting ....")

  methLevelOut = runParallel(methLevelMat, reshapeMat, binSize, covCut, runThread)

  outputFile = outputPrefix + ".txt"
  methLevelOut.to_csv(outputFile, sep = '\t', encoding='utf-8', index=False)

#####################
#### subfunction ####
#####################

def versions():
  verStr = "Program:\tbinMethLevel.py\nVersion:\t1.0"
  print(verStr)

def usageInfo():
  versions()
  print(__doc__)

def welcomeWords():
  welWords = "Welcome to use binMethLevel.py"
  print("*"*(len(welWords)+20))
  print("*"," "*(7),welWords," "*(7),"*")
  print("*"," "*(len(welWords)+16),"*")
  print("*"," "*(12),"Fei Sang,DeepSeq,QMC"," "*(12),"*")
  print("*"," "*(5),"Life Sciences,Univ. of Nottingham"," "*(6),"*")
  print("*"," "*(10),"Nottingham, NG7 2UH, UK"," "*(11),"*")
  print("*"," "*(9),"fei.sang@nottingham.ac.uk"," "*(10),"*")
  print("*"*(len(welWords)+20))

def runTime(strTmp):
  runningTime = str(datetime.now())
  runningTime = runningTime.split(".")[0]
  print(runningTime,">>>>",strTmp)

def reshapeMat(methLevel, binSize, covCut):
  bindex = 0
  chrID = ""
  methLevel_freq = np.array([])
  methLevel_covg = np.array([])
  methLevel_bin = pd.DataFrame()

  # reshape matrix
  pbar = tqdm(range(len(methLevel.index)), ncols = 100)
  for rindex in pbar:
    if chrID == methLevel.iloc[rindex,0]:
      if (methLevel.iloc[rindex,1] >= bindex) & (methLevel.iloc[rindex,1] < (bindex+binSize)):
        methLevel_freq = np.append(methLevel_freq, methLevel.iloc[rindex,4])
        methLevel_covg = np.append(methLevel_covg, methLevel.iloc[rindex,5])
      else:
        methLevel_bin = methLevel_bin.append({'chrbin'              : chrID+"_"+str(bindex)+"_"+str(bindex+binSize),
                                              'freq'                : ','.join(map(str,methLevel_freq)),
                                              'covg'                : ','.join(map(str,methLevel_covg)),
                                              'avgfreq'             : round(methLevel_freq.mean(),3),
                                              'numcov'+str(covCut)  : sum(methLevel_covg>=covCut),
                                              'numcov'              : np.size(methLevel_covg),
                                              'covfrac'             : round(sum(methLevel_covg>=covCut)/np.size(methLevel_covg),3)}, ignore_index=True)
        bindex = bindex + binSize
        methLevel_freq = np.array([])
        methLevel_covg = np.array([])
        methLevel_freq = np.append(methLevel_freq, methLevel.iloc[rindex,4])
        methLevel_covg = np.append(methLevel_covg, methLevel.iloc[rindex,5])
    else:
      if rindex != 0 :
        methLevel_bin = methLevel_bin.append({'chrbin'              : chrID+"_"+str(bindex)+"_"+str(bindex+binSize),
                                              'freq'                : ','.join(map(str,methLevel_freq)),
                                              'covg'                : ','.join(map(str,methLevel_covg)),
                                              'avgfreq'             : round(methLevel_freq.mean(),3),
                                              'numcov'+str(covCut)  : sum(methLevel_covg>=covCut),
                                              'numcov'              : np.size(methLevel_covg),
                                              'covfrac'             : round(sum(methLevel_covg>=covCut)/np.size(methLevel_covg),3)}, ignore_index=True)
      bindex = int(methLevel.iloc[rindex,1]/binSize)*binSize
      chrID = methLevel.iloc[rindex,0]
      methLevel_freq = np.array([])
      methLevel_covg = np.array([])
      methLevel_freq = np.append(methLevel_freq, methLevel.iloc[rindex,4])
      methLevel_covg = np.append(methLevel_covg, methLevel.iloc[rindex,5])
  pbar.close()

  # last element in the data frame
  methLevel_bin = methLevel_bin.append({'chrbin'              : chrID+"_"+str(bindex)+"_"+str(bindex+binSize),
                                        'freq'                : ','.join(map(str,methLevel_freq)),
                                        'covg'                : ','.join(map(str,methLevel_covg)),
                                        'avgfreq'             : round(methLevel_freq.mean(),3),
                                        'numcov'+str(covCut)  : sum(methLevel_covg>=covCut),
                                        'numcov'              : np.size(methLevel_covg),
                                        'covfrac'             : round(sum(methLevel_covg>=covCut)/np.size(methLevel_covg),3)}, ignore_index=True)
                    
  methLevel_bin = methLevel_bin[['chrbin','freq','covg','avgfreq','numcov'+str(covCut),'numcov','covfrac']]
  methLevel_bin[['numcov']] = methLevel_bin[['numcov']].astype(int)
  methLevel_bin[['numcov'+str(covCut)]] = methLevel_bin[['numcov'+str(covCut)]].astype(int)

  return methLevel_bin

def runParallel(df, func, b, c, t):
  df_split = np.array_split(df, t)
  pool = mp.Pool(t)
  df_out = pool.map(partial(func,binSize=b,covCut=c), df_split)
  pool.close()
  pool.join()
  return pd.concat(df_out)

#####################
#### program run ####
#####################

if __name__ == "__main__":
  main(sys.argv[1:])

