import sys
import getopt
from subprocess import call
from classes.MirnaDB import MiRNADB

# Parameters:
#
# -s: Representative species name
# -o: File name for positive feature set
opts, extraparams = getopt.getopt(sys.argv[1:], 's:o:')
for o,p in opts:
	if o == '-s':
		species = p
	if o == '-o':
		outPath = p

# Load up the mirbase database and clusters, export clustered features.
mirdb = MiRNADB()
mirdb.load_from_file('mirbase_19.mirdb')
mirdb.filter_by_cluster('mirbase_19.clus', inSpecies=species)
mirdb.filter_top_results(691)
mirdb.export_features('data/'+outPath)
mirdb.export_fasta('data/'+outPath+'.fasta')
mirdb.export_labels('data/'+outPath+'.labels')