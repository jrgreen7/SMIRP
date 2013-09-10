import getopt, sys
from classes.FeatureSet import FeatureSet

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:f:')
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-f':
		outFormat = p

inFormat = inPath.split('.')[-1]
noFormatName = ""
for text in inPath.split('.')[:-1]:
	noFormatName += text
	noFormatName += '.'
outPath = noFormatName + outFormat

if inFormat in ['micropred', 'features', 'huntmi', 'csv', 'svm', 'libsvm', 'arff']:
	fs = FeatureSet()
	fs.load(inPath)
	fs.export(outPath)