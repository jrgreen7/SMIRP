import os
import sys
import getopt
from classes.FeatureSet import FeatureSet
from subprocess import call

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		inPath = p


# Make fold data
call('progs/HuntMi/progs/RNAfold/RNAfold -noPS < data/'+inPath+' > data/'+inPath+'.fold', shell=True)

# Make dustmasker feature. Delete temporary files associated with dustmasker feature.
call('python progs/HuntMi/progs/additional_features/dustmasker.py data/'+inPath+' data/'+inPath+'.dustmasker', shell=True)
call('rm tmp', shell=True)
call('rm data/'+inPath+'.dustmasker.tmp', shell=True)

# Make loops, translate, and triplet feature files
call('python progs/HuntMi/progs/additional_features/loops.py data/'+inPath+'.fold data/'+inPath+'.loops', shell=True)
call('python progs/HuntMi/progs/additional_features/translate.py data/'+inPath+' data/'+inPath+'.translate', shell=True)
call('python progs/HuntMi/progs/additional_features/triplet.py data/'+inPath+'.fold data/'+inPath+'.triplet', shell=True)

fs = FeatureSet()
fs.load_micropred('data/'+inPath+'.micropred')
fs.add_features_from_micropred('data/'+inPath+".dustmasker")
fs.add_features_from_micropred('data/'+inPath+'.triplet')
fs.add_features_from_micropred('data/'+inPath+'.loops')
fs.add_features_from_micropred('data/'+inPath+'.translate')
fs.export_micropred('data/'+inPath+'.huntmi')

call('rm data/'+inPath+'.fold', shell=True)
call('rm data/'+inPath+'.translate', shell=True)
call('rm data/'+inPath+'.triplet', shell=True)
call('rm data/'+inPath+'.loops', shell=True)
call('rm data/'+inPath+'.dustmasker', shell=True)