import sys
import getopt
from classes.ResultSet import ResultSet

# Parameters:
#
# -f: Fasta file containing hairpin sequences
# -r: Micropred result file
# -o: output file
opts, extraparams = getopt.getopt(sys.argv[1:], 'f:r:o:')
for o,p in opts:
	if o == '-f':
		fastaPath = p
	if o == '-r':
		resultPath = p
	if o == '-o':
		outPath = p

rs = ResultSet()
rs.load('data/'+fastaPath, 'data/'+resultPath)
rs.sort()
rs.export('data/'+outPath)