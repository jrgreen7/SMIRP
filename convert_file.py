import getopt, sys
from classes.FeatureSet import FeatureSet

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:f:')
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-f':
		outFormat = p

inFormat = inPath.rsplit('.',1)[0]
noFormatName = inPath.rsplit('.',1)[1]:

outPath = noFormatName + outFormat

if inFormat in ['micropred', 'features', 'huntmi', 'csv', 'svm', 'libsvm', 'arff']:
	fs = FeatureSet()
	fs.load(inPath)
	fs.export(outPath)