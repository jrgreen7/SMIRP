from subprocess import call
import sys
import getopt
import classes.FastaOperations as FastaOps
import classes.FoldOperations as FoldOps
import classes.FileConversion as FileConversion
from classes.SequenceList import *

# Parameters:
# -i: Input fasta file
#
# Output is 
opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		inPath = p

print "Re-formatting fasta."
# Step one: turn the fasta into something that RNALfold will work with
FastaOps.remove_newlines('data/'+inPath, 'data/'+inPath+'.fixed')

print "Finding all folds with RNALfold."
call('progs/ViennaRNA-1.8.5/Progs/RNALfold -d2 -noLP -L 120 < data/'+inPath+'.fixed > data/'+inPath+'.folds', shell=True)

print "Filtering folds down to hairpins."
FoldOps.filter_hairpins('data/'+inPath+'.folds', 'data/'+inPath+'.hairpins')

FileConversion.RNAL_to_fasta('data/'+inPath+'.hairpins', 'data/folds_from_'+inPath)

print "Removing redundant hairpins."
sl = SequenceList()
sl.load_fasta('data/folds_from_'+inPath)
sl.remove_all_redundant()
sl.export_fasta('data/'+inPath+'.nr.hairpins')
print "Done!"