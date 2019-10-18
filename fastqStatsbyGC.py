#################################################################
#								#
#  This script is for statistical analysis for reads in fasta.	#
#								#
#################################################################

#!/usr/bin/python

##########
# ReadMe #
##########

'''
USAGE:\tfastqStatsbyGC.py -g -p -i <inputfile> -o <outputfile prefix>

Arguments:
\t-g, --gc              GC cutoff (default: 40%)
\t-p, --paired          pair ended reads
\t-i, --input           input file in fastq format
\t-o, --outputprefix    output file prefix

Flags:
\t-h, --help            help information
\t-v, --version         version information

Example:
python fastqStatsbyGC.py -i singles.fastq -o test
python fastqStatsbyGC.py -g 40 -p -i pairs.read_1.fastq,pairs.read_2.fastq -o test

Developer:
\tFei Sang
\tDeep Seq, Queen's Medical School
\tSchool of Life Sciences
\tUniversity of Nottingham
\tEmail: fei.sang@nottingham.ac.uk
''' 

#################
# main function #
#################

import sys, getopt
from datetime import datetime

import math

from Bio import SeqIO
from Bio.SeqUtils import GC

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
from scipy.stats import gaussian_kde
from itertools import izip

def main(argvs):

	# parameters
	gcCut = 40
	peRead = False
	inputFile = ''
	outputPrefix = ''

	try:
		opts, args = getopt.getopt(argvs,"vhg:pi:o:",["version","help","gc","paired","input=","outputprefix="])
		if len(opts) == 0:
			print "Please use the correct arguments, for usage type -h"
			sys.exit(2)
	except getopt.GetoptError, err:
		print str(err)
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-v", "--version"):
			verbose = True
			versions()
			sys.exit()
		elif opt in ("-h", "--help"):
			usageInfo()
			sys.exit()
		elif opt in ("-g", "--gc"):
			gcCut = arg
		elif opt in ("-p", "--paired"):
			peRead = True
		elif opt in ("-i", "--input"):
			inputFile = arg
		elif opt in ("-o", "--outputprefix"):
			outputPrefix = arg
		else:
			assert False, "unhandled option"

	if gcCut == "": gcCut = 40
	if inputFile == "": sys.exit("Please use the correct arguments, option -i missing")
	if outputPrefix == "": sys.exit("Please use the correct arguments, option -o missing")

	# processing
	welcomeWords()

	runningTime = str(datetime.now())
	runningTime = runningTime.split(".")[0]
	print runningTime,">>>>","Fetch read information & build index"

	processIndex = 0
	processPrint = ""

	if peRead:
		firstReadFile = inputFile.split(",")[0]
		secondReadFile = inputFile.split(",")[1]

		firstR1OutFile = outputPrefix + '.GCabove' + str(gcCut) + '.read_1.fastq'
		secondR1OutFile = outputPrefix + '.GCbelow' + str(gcCut) + '.read_1.fastq'
		firstR2OutFile = outputPrefix + '.GCabove' + str(gcCut) + '.read_2.fastq'
		secondR2OutFile = outputPrefix + '.GCbelow' + str(gcCut) + '.read_2.fastq'

		firstR1OF = open(firstR1OutFile, "w")
		secondR1OF = open(secondR1OutFile, "w")
		firstR2OF = open(firstR2OutFile, "w")
		secondR2OF = open(secondR2OutFile, "w")

		secondReadFileDict = SeqIO.index(secondReadFile, "fastq")

		for mR1Info in SeqIO.parse(firstReadFile, "fastq"):
			mR2Info = mR1Info.id
			if (GC(mR1Info.seq) > int(gcCut)) and (GC(secondReadFileDict[mR2Info].seq) > int(gcCut)): 
				SeqIO.write(mR1Info, firstR1OF, "fastq")
				SeqIO.write(secondReadFileDict[mR2Info], firstR2OF, "fastq")
			elif (GC(mR1Info.seq) <= int(gcCut)) and (GC(secondReadFileDict[mR2Info].seq) <= int(gcCut)):
				SeqIO.write(mR1Info, secondR1OF, "fastq")
				SeqIO.write(secondReadFileDict[mR2Info], secondR2OF, "fastq")

			processIndex = processIndex + 1
			processPrint = runningTime + ' ' + ">>>> Processing:" + ' ' + str(processIndex) + '\r'
			if processIndex%1000 == 0:
				sys.stdout.write(processPrint)
				sys.stdout.flush()

		print "\n"

		outputPrefixTmp = outputPrefix + '.read_1'
		statPlot(firstReadFile, outputPrefixTmp, gcCut)
		outputPrefixTmp = outputPrefix + '.read_2'
		statPlot(secondReadFile, outputPrefixTmp, gcCut)

		secondR1OF.close()
		firstR2OF.close()
		secondR2OF.close()
		firstR1OF.close()

	else:
		firstOutFile = outputPrefix + '.GCabove' + str(gcCut) + '.fastq'
		secondOutFile = outputPrefix + '.GCbelow' + str(gcCut) + '.fastq'

		firstOF = open(firstOutFile, "w")
		secondOF = open(secondOutFile, "w")

		for mReadInfo in SeqIO.parse(inputFile, "fastq"):
			if GC(mReadInfo.seq) > int(gcCut): 
				SeqIO.write(mReadInfo, firstOF, "fastq")
			else:
				SeqIO.write(mReadInfo, secondOF, "fastq")

			processIndex = processIndex + 1
			processPrint = runningTime + ' ' + ">>>> Processing:" + ' ' + str(processIndex) + '\r'
			if processIndex%1000 == 0:
				sys.stdout.write(processPrint)
				sys.stdout.flush()

		print "\n"

		statPlot(inputFile, outputPrefix, gcCut)

		secondOF.close()
		firstOF.close()

###############
# subfunction #
###############

def versions():
	verStr = "Program:\tfastaStats.py\nVersion:\t1.0"
	print verStr

def usageInfo():
	versions()
	print (__doc__)

def welcomeWords():
	welWords = "Welcome to use fastaStats.py"
	print "*"*(len(welWords)+20)
	print "*"," "*(7),welWords," "*(7),"*"
	print "*"," "*(len(welWords)+16),"*"
	print "*"," "*(11),"Fei Sang,DeepSeq,QMC"," "*(11),"*"
	print "*"," "*(4),"Life Sciences,Univ. of Nottingham"," "*(5),"*"
	print "*"," "*(9),"Nottingham, NG7 2UH, UK"," "*(10),"*"
	print "*"," "*(8),"fei.sang@nottingham.ac.uk"," "*(9),"*"
	print "*"*(len(welWords)+20)

def statPlot(subInput, subOutput, subGC):
	readLength = list()
	readGC = list()

	readLengthGroup1 = list()
	readGCGroup1 = list()
	readLengthGroup2 = list()
	readGCGroup2 = list()

	subProcessIndex = 0
	subProcessPrint = ""
	
	subRunningTime = str(datetime.now())
	subRunningTime = subRunningTime.split(".")[0]
	print subRunningTime,">>>>","Fetch read information for plotting"

	for readInfo in SeqIO.parse(subInput, "fastq"):
		readLength.append(math.log10(len(readInfo.seq)))
		readGC.append(GC(readInfo.seq))
		if GC(readInfo.seq) > int(subGC): 
			readLengthGroup1.append(math.log10(len(readInfo.seq)))
			readGCGroup1.append(GC(readInfo.seq))
		else:
			readLengthGroup2.append(math.log10(len(readInfo.seq)))
			readGCGroup2.append(GC(readInfo.seq))

		subProcessIndex = subProcessIndex + 1
		subProcessPrint = subRunningTime + ' ' + ">>>> Processing:" + ' ' + str(subProcessIndex) + '\r'
		if subProcessIndex%1000 == 0:
			sys.stdout.write(subProcessPrint)
			sys.stdout.flush()

	print "\n"

	subRunningTime = str(datetime.now())
	subRunningTime = subRunningTime.split(".")[0]
	print subRunningTime,">>>>","Plot distribution of read length, GC content"

	#############
	# all reads #
	#############
	readLengthDensity = gaussian_kde(readLength)
	readGCDensity = gaussian_kde(readGC)

	readLengthXS = numpy.linspace(min(readLength),max(readLength),100)
	readGCXS = numpy.linspace(min(readGC),max(readGC),100)

	pltFig = plt.figure()
	pltObjA = pltFig.add_subplot(1,1,1)
	pltObjA.hist(readLength, bins=100, color='lightgrey')
	pltObjB = pltObjA.twinx()
	pltObjB.plot(readLengthXS, readLengthDensity(readLengthXS), color='red', linewidth=2)
	pltObjA.set_xlabel("Contig Length (log)")
	pltObjA.set_ylabel("Frequency")
	pltObjB.set_ylabel("Density")
	titleFig = 'Length Distribution (All Reads)'
	plt.title(titleFig)
	nameFig = subOutput + '.All.readLengthDistri.jpg'
	pltFig.savefig(nameFig)

	pltFig = plt.figure()
	pltObjA = pltFig.add_subplot(1,1,1)
	pltObjA.hist(readGC, bins=100, color='lightgrey')
	pltObjB = pltObjA.twinx()
	pltObjB.plot(readGCXS, readGCDensity(readGCXS), color='red', linewidth=2)
	pltObjA.set_xlabel("Contig GC Content (%)")
	pltObjA.set_ylabel("Frequency")
	pltObjB.set_ylabel("Density")
	titleFig = 'GC Content Distribution (All Reads)'
	plt.title(titleFig)
	nameFig = subOutput + '.All.readGCDistri.jpg'
	pltFig.savefig(nameFig)

	################
	# group1 reads #
	################
	readLengthDensity = gaussian_kde(readLengthGroup1)
	readGCDensity = gaussian_kde(readGCGroup1)

	readLengthXS = numpy.linspace(min(readLengthGroup1),max(readLengthGroup1),100)
	readGCXS = numpy.linspace(min(readGCGroup1),max(readGCGroup1),100)

	pltFig = plt.figure()
	pltObjA = pltFig.add_subplot(1,1,1)
	pltObjA.hist(readLengthGroup1, bins=100, color='lightgrey')
	pltObjB = pltObjA.twinx()
	pltObjB.plot(readLengthXS, readLengthDensity(readLengthXS), color='red', linewidth=2)
	pltObjA.set_xlabel("Contig Length (log)")
	pltObjA.set_ylabel("Frequency")
	pltObjB.set_ylabel("Density")
	titleFig = 'Length Distribution' + ' (' + 'GC above ' + str(subGC) + '%)'
	plt.title(titleFig)
	nameFig = subOutput + '.GCabove' + str(subGC) + '.readLengthDistri.jpg'
	pltFig.savefig(nameFig)

	pltFig = plt.figure()
	pltObjA = pltFig.add_subplot(1,1,1)
	pltObjA.hist(readGCGroup1, bins=100, color='lightgrey')
	pltObjB = pltObjA.twinx()
	pltObjB.plot(readGCXS, readGCDensity(readGCXS), color='red', linewidth=2)
	pltObjA.set_xlabel("Contig GC Content (%)")
	pltObjA.set_ylabel("Frequency")
	pltObjB.set_ylabel("Density")
	titleFig = 'GC Content Distribution' + ' (' + 'GC above ' + str(subGC) + '%)'
	plt.title(titleFig)
	nameFig = subOutput + '.GCabove' + str(subGC) + '.readGCDistri.jpg'
	pltFig.savefig(nameFig)

	################
	# group2 reads #
	################
	readLengthDensity = gaussian_kde(readLengthGroup2)
	readGCDensity = gaussian_kde(readGCGroup2)

	readLengthXS = numpy.linspace(min(readLengthGroup2),max(readLengthGroup2),100)
	readGCXS = numpy.linspace(min(readGCGroup2),max(readGCGroup2),100)

	pltFig = plt.figure()
	pltObjA = pltFig.add_subplot(1,1,1)
	pltObjA.hist(readLengthGroup2, bins=100, color='lightgrey')
	pltObjB = pltObjA.twinx()
	pltObjB.plot(readLengthXS, readLengthDensity(readLengthXS), color='red', linewidth=2)
	pltObjA.set_xlabel("Contig Length (log)")
	pltObjA.set_ylabel("Frequency")
	pltObjB.set_ylabel("Density")
	titleFig = 'Length Distribution' + ' (' + 'GC below ' + str(subGC) + '%)'
	plt.title(titleFig)
	nameFig = subOutput + '.GCbelow' + str(subGC) + '.readLengthDistri.jpg'
	pltFig.savefig(nameFig)

	pltFig = plt.figure()
	pltObjA = pltFig.add_subplot(1,1,1)
	pltObjA.hist(readGCGroup2, bins=100, color='lightgrey')
	pltObjB = pltObjA.twinx()
	pltObjB.plot(readGCXS, readGCDensity(readGCXS), color='red', linewidth=2)
	pltObjA.set_xlabel("Contig GC Content (%)")
	pltObjA.set_ylabel("Frequency")
	pltObjB.set_ylabel("Density")
	titleFig = 'GC Content Distribution' + ' (' + 'GC below ' + str(subGC) + '%)'
	plt.title(titleFig)
	nameFig = subOutput + '.GCbelow' + str(subGC) + '.readGCDistri.jpg'
	pltFig.savefig(nameFig)

###############
# program run #
###############

if __name__ == "__main__":
	main(sys.argv[1:])

