import sys
import getopt

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		inPath = p

with open('data/selected.'+inPath+'.-21.features', 'r') as inFile:
	with open('data/'+inPath+'.micropred', 'w') as outFeatures:
		with open('data/'+inPath+'.labels', 'w') as outLabels:
			for line in inFile:
				if line[0] == '>':
					outLabels.write(line)
				if line[0] in '1234567890-.':
					outFeatures.write(line)
