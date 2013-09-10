import sys
import getopt
from subprocess import call
from classes.FeatureSet import FeatureSet

# Parameters:
#
# -p: File name for positive feature set (any file type)
# -n: File name for negative feature set (any file type)
# -o: Name of output LibSVM model

opts, extraparams = getopt.getopt(sys.argv[1:], 'o:p:n:')
for o,p in opts:
	if o == '-p':
		posPath = p
	if o == '-n':
		negPath = p
	if o == '-o':
		outPath = p

# Aggregate inputs, export to libsvm file
fs = FeatureSet()
fs.load('data/'+posPath, patternClass = 'real')
fs.add_instances('data/'+negPath, patternClass = 'pseudo')
fs.weka_smote()
fs.libsvm_scale(paramOut = 'models/'+outPath+'.scale')
fs.export('tmp.libsvm')
# Build model
call('svm-train -c 10000000 -d 1 -h 1 -e 0.001 -g 0.0019531 -b 1 tmp.libsvm models/'+outPath+'.model', shell=True)
# Clean up
call('rm tmp.libsvm', shell=True)