import classes.FastaOperations as fo
from classes.SequenceList import *

# fo.split_fasta('data/folds_from_AHGY01.fa', 10)
for i in range(10):
	sl = SequenceList()
	sl.load_fasta('data/folds_from_AHGY01.fa.'+str(i))
	sl.remove_redundant()
	sl.export_fasta('data/AHGY01.fa.nr.hairpins'+'.'+str(i))
	print str((i+1)*10)+'% complete removing redundant hairpins'
fo.merge_fasta('data/AHGY01.fa.nr.hairpins', 10)