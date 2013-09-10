import sys
import getopt
from subprocess import call
from classes.FeatureSet import FeatureSet
from classes.SequenceList import SequenceList

# Parameters:
#
# -s: Representative species name
# -c: Coding region filename
# -i: Input genome
# -o: File name for positive feature set
# -n number of cores handling workload
opts, extraparams = getopt.getopt(sys.argv[1:], 's:c:i:o:n:')
for o,p in opts:
	if o == '-s':
		species = p
	if o == '-c':
		negPath = p
	if o == '-i':
		inPath = p
	if o == '-o':
		outPath = p
	if o == '-n':
		numThreads = p

################################################
# Build positive training set
################################################
print "### Building positive set from miRbase database"
speciesFilename = species.replace(' ', '_')
call('python build_positiveset.py -s "'+species+'" -o '+speciesFilename+'.features', shell=True)

################################################
# Build negative training set
################################################
print "### Building negative set hairpins from coding regions"
call('python extract_hairpins.py -i '+negPath, shell=True)
print "### Extracting micropred features from coding regions"
sl = SequenceList()
sl.load_fasta('data/'+negPath+'.nr.hairpins')
sl.select_random(10000)
sl.export_fasta('data/'+negPath+'.nr.hairpins')
call('python build_micropred_features.py -i '+negPath+'.nr.hairpins -n '+numThreads, shell=True)
# call('python build_huntmi_features.py -i '+negPath+'.nr.hairpins')

################################################
# Build LibSVM model
################################################
print "### Building LibSVM model"
call('python build_model.py -p '+speciesFilename+'.features -n '+negPath+'.nr.hairpins.micropred -o '+speciesFilename, shell=True)

################################################
# Build feature set from hairpin candidates in genome of interest
################################################
print "### Building hairpins from genome under exploration"
call('python extract_hairpins.py -i '+inPath, shell=True)
print "### Extracting micropred features from genome under exploration"
call('python build_micropred_features.py -i '+inPath+'.nr.hairpins -n '+numThreads, shell=True)

################################################
# Run svm-predict on all hairpin candidates in genome of interest
################################################
fs = FeatureSet()
fs.load('data/'+inPath+'.nr.hairpins.micropred', patternClass = 'real')
fs.libsvm_scale(params='models/'+speciesFilename+'.scale')
fs.export('data/'+inPath+'.nr.hairpins.libsvm')
call('progs/libsvm-3.14/svm-predict -b 1 data/'+inPath+'.nr.hairpins.libsvm models/'+speciesFilename+'.model data/'+inPath+'.nr.hairpins.results', shell=True)