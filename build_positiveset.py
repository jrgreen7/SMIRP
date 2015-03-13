import sys
import getopt
from subprocess import call
from classes.MirnaDB import MiRNADB

# Parameters:
#
# -s: Representative species name
# -o: File name for positive feature set
# -n: Number of miRNA in positive sset
opts, extraparams = getopt.getopt(sys.argv[1:], 's:o:n:')
num = 691
for o,p in opts:
	if o == '-s':
		species = p
	if o == '-o':
		outPath = p
	if o == '-n':
		num = int(p)

# Load up the mirbase database and clusters, export clustered features.
mirdb = MiRNADB()
mirdb.load_from_file('mirbase_21.mirdb')
mirdb.filter_by_cluster('mirbase_21.clus', inSpecies=species)
mirdb.filter_top_results(num)
mirdb.export_features('data/'+outPath)
mirdb.export_fasta('data/'+outPath+'.fasta')
# mirdb.export_labels('data/'+outPath+'.labels')