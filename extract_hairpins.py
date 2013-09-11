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

# Parameters:
# -i: Input fasta file
# -n: Number of threads for multi-threaded use
#
# Output is a non-redundant hairpins file, <inpath>.nr.hairpins
opts, extraparams = getopt.getopt(sys.argv[1:], 'i:n:')
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-n':
		numThreads = int(p)

class myThread(threading.Thread):
	def __init__(self, inPath):
		threading.Thread.__init__(self)
		self.inPath = inPath

	def run(self):
		print "Finding all folds in "+self.inPath+" with RNALfold."
		call('progs/ViennaRNA-1.8.5/Progs/RNALfold -d2 -noLP -L 120 < data/tmp/'+self.inPath+' > data/tmp/'+self.inPath+'.folds', shell=True)
		print "Filtering folds  in "+self.inPath+" down to hairpins."
		FoldOps.filter_hairpins('data/tmp/'+self.inPath+'.folds', 'data/tmp/'+self.inPath+'.hairpins')
		FileConversion.RNAL_to_fasta('data/tmp/'+self.inPath+'.hairpins', 'data/tmp/folds_from_'+self.inPath)
		print "Removing redundant hairpins from "+self.inPath+"."
		sl = SequenceList()
		sl.load_fasta('data/tmp/folds_from_'+self.inPath)
		sl.remove_all_redundant()
		sl.export_fasta('data/tmp/'+self.inPath+'nrhairpins')


print "Re-formatting fasta."
# Step one: turn the fasta into something that RNALfold will work with
FastaOps.remove_newlines('data/'+inPath, 'data/tmp/'+inPath+'.fixed')
# Step two: split the fasta for mutli-threaded processing
FastaOps.split_fasta('data/tmp/'+inPath+'.fixed', numThreads)
# Step three: Launch threads
threadPath = inPath
threadExt = 'fixed'
threads = []
for i in range(numThreads):
	threads.append(myThread(threadPath+'.'+str(i)+'.'+threadExt))
for thread in threads:
	thread.start()
for thread in threads:
	thread.join()

FastaOps.merge_fasta('data/tmp/'+inPath+'.fixednrhairpins', numThreads)
call('cp data/tmp/'+inPath+'.fixednrhairpins data/'+inPath+'.nr.hairpins', shell=True)
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