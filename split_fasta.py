import sys
import getopt
import classes.FastaOperations as fo

# Arguments
############
# -i: Input file (in data/ directory)
# -n: Number of output files

# Output
#########
# n fasta files, with 


opts, extraparams = getopt.getopt(sys.argv[1:], 'i:n:')
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-n':
		numFiles = p

fo.split_fasta('data/'+inPath, int(numFiles))