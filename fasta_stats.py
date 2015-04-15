import sys
import subprocess
import getopt
from classes.ResultSet import ResultSet
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm

# Parameters:
# -f: Input fasta file containing sequences
# -h: HMP (not hmp20) file 
opts, extraparams = getopt.getopt(sys.argv[1:], 'f:h:')
for o,p in opts:
	if o == '-f':
		inPath = p
	if o == '-h':
		featPath = p


subprocess.call('mkdir results/'+inPath, shell=True)

headers = [x.strip() for x in open('data/'+inPath, 'r').readlines()[0::2]]
seqs = [x.strip() for x in open('data/'+inPath, 'r').readlines()[1::2]]
MFEs = []
with open('data/'+featPath, 'r') as inFeats:
	for line in inFeats:
		MFEs.append(float(line.split(',')[0]))

def get28GCs(seqs):
	GCs = []
	fasta = open('data/'+inPath, 'r').readlines()
	for line in seqs:
		seq = line.strip()[1:8]
		seqLen = float(len(seq))
		seqGC = (float(seq.count('G')) + float(seq.count('C'))) / seqLen
		GCs.append(seqGC)
	return GCs

def getGCs(seqs):
	GCs = []
	fasta = open('data/'+inPath, 'r').readlines()
	for line in seqs:
		seq = line.strip()
		seqLen = float(len(seq))
		seqGC = (float(seq.count('G')) + float(seq.count('C'))) / seqLen
		GCs.append(seqGC)
	return GCs

def getLens(seqs):
	lens = []
	fasta = open('data/'+inPath, 'r').readlines()
	for line in seqs:
		seq = line.strip()
		seqLen = float(len(seq))
		lens.append(seqLen)
	return lens

def getGCratios(seqs):
	ratios = []
	fasta = open('data/'+inPath, 'r').readlines()
	for line in seqs:
		seq = line.strip()
		seqLen = float(len(seq))
		seqGC = (float(seq.count('G')) + float(seq.count('C'))) / seqLen
		seqG = float(seq.count('G')) / seqLen
		ratios.append(seqG / seqGC)
	return ratios


def getAUratios(seqs):
	ratios = []
	fasta = open('data/'+inPath, 'r').readlines()
	for line in seqs:
		seq = line.strip()
		seqLen = float(len(seq))
		seqAU = (float(seq.count('A')) + float(seq.count('U'))) / seqLen
		seqA = float(seq.count('A')) / seqLen
		ratios.append(seqA / seqAU)
	return ratios

GCs = getGCs(seqs)
GC28s = get28GCs(seqs)
lens = getLens(seqs)
AAUs = getAUratios(seqs)
GGCs = getGCratios(seqs)

outFile = open("results/"+inPath+"/stats.csv", 'w')

outFile.write("Header,Sequence,Length,\"GC content\",\"GC content (2-8)\",G/GC,A/AU\n")

for idx, seq in enumerate(seqs):
	outFile.write("\""+headers[idx]+"\""+',')
	outFile.write(seq+',')
	outFile.write(str(lens[idx])+',')
	outFile.write(str(GCs[idx])+',')
	outFile.write(str(GC28s[idx])+',')
	outFile.write(str(GGCs[idx])+',')
	outFile.write(str(AAUs[idx])+'\n')
outFile.write("Header,Sequence,Length,\"GC content\",\"GC content (2-8)\",G/GC,A/AU\n")
outFile.write(",Average,")
outFile.write(str(sum(lens)/len(lens))+',')
outFile.write(str(sum(GCs)/len(GCs))+',')
outFile.write(str(sum(GC28s)/len(GC28s))+',')
outFile.write(str(sum(GGCs)/len(GGCs))+',')
outFile.write(str(sum(AAUs)/len(AAUs))+'\n')
outFile.close()

# GCmean = sum(GCs)/len(GCs)
# otherGCs = getGCs("Physarum_miRNA.seqs")
# otherMean = sum(otherGCs)/len(otherGCs)
# print GCmean

# (mu, sigma) = norm.fit(GCs)
# print "Mu:", mu
# print "Sigma:", sigma

# n, bins, patches = plt.hist(GCs, normed=1, bins=50, color='c')
# # plt.hist(otherGCs, normed=1, bins=5, color='r', alpha = 0.5)
# y = mlab.normpdf(bins, mu, sigma)
# l = plt.plot(bins, y, 'b--', linewidth=2)
# plt.axvline(GCmean, color='b', linestyle='dashed', linewidth=2)
# for pt in otherGCs:
# 	plt.axvline(pt, color='g', linewidth=1)
# plt.show()

