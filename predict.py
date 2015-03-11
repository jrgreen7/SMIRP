import sys
import getopt
from subprocess import call
from classes.FeatureSet import FeatureSet
from classes.ResultSet import ResultSet

# Parameters:
#
# -m: model name (<-m>.scale and <-m>.model should exist in models directory)
# -i: File containing input feature data
opts, extraparams = getopt.getopt(sys.argv[1:], 'm:i:')
for o,p in opts:
	if o == '-m':
		modelName = p
	if o == '-i':
		inPath = p

fs = FeatureSet()
fs.load('data/'+inPath, patternClass = 'real')
fs.libsvm_scale(params='models/'+modelName+'.scale')
fs.export('tmp.libsvm')
call('progs/libsvm-3.14/svm-predict -b 1 tmp.libsvm models/'+modelName+'.model data/'+inPath+'.results', shell=True)
# call('rm tmp.libsvm', shell=True)