from subprocess import call
import threading
import sys
import getopt
import os.path
import os
import classes.FastaOperations as FastaOps
import classes.FoldOperations as FoldOps
import classes.FileConversion as FileConversion
from classes.SequenceList import *
import random

# Parameters:
# -i: Input fasta file
# -n: Number of threads for multi-threaded use
#
# Output is a non-redundant hairpins file, <inpath>.nr.hairpins
opts, extraparams = getopt.getopt(sys.argv[1:], 'i:n:l:m:p:h:')
hairpinLength = 100
basePairs = 18
minMFE = -15.00
numHairpins = 0
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-n':
		numThreads = int(p)
	if o == '-l':
		hairpinLength = int(p)
	if o == '-m':
		minMFE = float(p)
	if o == '-p':
		basePairs = int(p)
	if o == '-h':
		numHairpins = int(p)



class myThread(threading.Thread):
	def __init__(self, inPath):
		threading.Thread.__init__(self)
		self.inPath = inPath

	def run(self):
		call('progs/ViennaRNA-2.1.8/Progs/RNALfold -d2 --noLP -L '+str(hairpinLength)+' < data/tmp/'+self.inPath+' > data/tmp/'+self.inPath+'.folds', shell=True)
		# call('progs/ViennaRNA-1.8.5/Progs/RNAfold < data/tmp/'+self.inPath+' > data/tmp/'+self.inPath+'.folds', shell=True)
		FoldOps.filter_hairpins('data/tmp/'+self.inPath+'.folds', 'data/tmp/'+self.inPath+'.hairpins', minMFE, basePairs)
		FileConversion.RNAL_to_fasta('data/tmp/'+self.inPath+'.hairpins', 'data/tmp/folds_from_'+self.inPath)
		sl = SequenceList()
		sl.load_fasta('data/tmp/folds_from_'+self.inPath)
		sl.remove_all_redundant()
		sl.export_fasta('data/tmp/'+self.inPath+'nrhairpins')


# Step one: turn the fasta into something that RNALfold will work with
FastaOps.remove_newlines('data/'+inPath, 'data/tmp/'+inPath+'.fixed')
FastaOps.convert_DNA_to_RNA('data/tmp/'+inPath+'.fixed', 'data/tmp/'+inPath+'.rna')
# Step two: split the fasta for mutli-threaded processing
FastaOps.split_fasta('data/tmp/'+inPath+'.rna', numThreads)
# Step three: Launch threads
threadPath = inPath
threadExt = 'rna'
threads = []
for i in range(numThreads):
	threads.append(myThread(threadPath+'.'+str(i)+'.'+threadExt))
for thread in threads:
	thread.start()
for thread in threads:
	thread.join()

FastaOps.merge_fasta('data/tmp/'+inPath+'.rnanrhairpins', numThreads)
FastaOps.remove_AU('data/tmp/'+inPath+'.rnanrhairpins', 'data/tmp/'+inPath+'.hairpins.noAU', 3)

call('cp data/tmp/'+inPath+'.hairpins.noAU data/'+inPath+'.nr.hairpins', shell=True)

if numHairpins != 0:
	outFile = open('data/'+inPath+'.nr.hairpins.'+str(numHairpins), 'w')
	inLines = open('data/'+inPath+'.nr.hairpins', 'r').readlines()
	inData = []
	for i in range(0,len(inLines)-2,2):
		inData.append(inLines[i]+inLines[i+1])

	if numHairpins < len(inData):
		outData = random.sample(inData, numHairpins)
	else:
		outData = inData
	for d in outData:
		outFile.write(d)




print "DONE!"

# print "Finding all folds in "+self.inPath+" with RNALfold."
# call('progs/ViennaRNA-1.8.5/Progs/RNALfold -d2 -noLP -L 120 < data/tmp/'+inPath+'.fixed > data/tmp/'+inPath+'.folds', shell=True)

# print "Filtering folds  in "+self.inPath+" down to hairpins."
# FoldOps.filter_hairpins('data/tmp/'+inPath+'.folds', 'data/tmp/'+inPath+'.hairpins')

# FileConversion.RNAL_to_fasta('data/tmp/'+inPath+'.hairpins', 'data/tmp/folds_from_'+inPath)

# print "Removing redundant hairpins from "+self.inPath+"."
# sl = SequenceList()
# sl.load_fasta('data/tmp/folds_from_'+inPath)
# sl.remove_all_redundant()
# sl.export_fasta('data/'+inPath+'.nr.hairpins')
# print "Done!"