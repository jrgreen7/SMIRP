import sys
import getopt
from classes.ResultSet import ResultSet

# Parameters:
#
# -f: Fasta file containing hairpin sequences
# -r: Micropred/HeteroMirPred result file
# -o: output file
opts, extraparams = getopt.getopt(sys.argv[1:], 'f:r:o:h:')
for o,p in opts:
	if o == '-f':
		fastaPath = p
	if o == '-r':
		resultPath = p
	if o == '-o':
		outPath = p
	if o == '-h':
		headerPath = p

rs = ResultSet()
rs.load_headers('data/'+fastaPath, 'data/'+resultPath, 'data/'+headerPath)
rs.sort()
rs.export('data/'+outPath)