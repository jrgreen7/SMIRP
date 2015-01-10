from subprocess import call
import sys
import getopt
import string
import random
import math
from collections import namedtuple
from operator import attrgetter
from classes.FeatureSet import FeatureSet

# Parameters:
# -p Filename for positive training set fasta file
# -n Filename for negative training set fasta file
# -f Number of folds used during outer-CV training
# -h prefix for hold-out data set (<h>_positive.libsvm and <-h>_negative.libsvm should exist in data/
opts, extraparams = getopt.getopt(sys.argv[1:], 'p:n:f:h:b:')
for o,p in opts:
	if o == '-p':
		posFile = p
	if o == '-n':
		negFile = p
	if o == '-f':
		numFolds = int(p)
	if o == '-h':
		hoSpec = p

def test_against_holdout(species, numFolds, posFile):
	outFile = open('data/'+species+'_roc_'+posFile.split('.')[0]+'.tsv', 'w')
	data = []
	Datum = namedtuple('Datum', 'predicted, true, probability')
	RocPoint = namedtuple('RocPoint', 'oneminusspec, sens')
	call('svm-scale -r data/params data/'+species+'_positive.libsvm > data/'+species+'_positive_scale.libsvm', shell=True)
	call('svm-scale -r data/params data/'+species+'_negative.libsvm > data/'+species+'_negative_scale.libsvm', shell=True)
	for i in range(numFolds):
		call('svm-predict -b 1 data/'+species+'_positive_scale.libsvm models/'+str(i)+'.model data/'+species+'_pos_'+str(i)+'.results', shell=True)
		call('svm-predict -b 1 data/'+species+'_negative_scale.libsvm models/'+str(i)+'.model data/'+species+'_neg_'+str(i)+'.results', shell=True)
	
	resultList = []
	for i in range(numFolds):
		truePos = 0.0
		trueNeg = 0.0
		falsePos = 0.0
		falseNeg = 0.0
		posFile = open('data/'+species+'_pos_'+str(i)+'.results', 'r')
		negFile = open('data/'+species+'_neg_'+str(i)+'.results', 'r')
		for line in posFile:
			if line[0] != 'l':
				data.append( Datum(predicted=line.split()[0], true='1', probability=float(line.split()[1])) )
				if line[0] == '1':
					truePos += 1.0
				else:
					falseNeg += 1.0
		for line in negFile:
			if line[0] != 'l':
				data.append( Datum(predicted=line.split()[0], true='-1', probability=float(line.split()[1])) )
				if line[0] == '1':
					falsePos += 1.0
				else:
					trueNeg += 1.0


		resultList.append( (truePos/(truePos+falseNeg),trueNeg/(trueNeg+falsePos)) )

	rocPoints = []
	data = sorted(data, key=attrgetter('probability'), reverse=True)
	with open('data/'+species+'_data_points.tsv', 'w') as outData:
		for point in data:
			outData.write(str(point.predicted)+'\t'+str(point.true)+'\t'+str(point.probability)+'\n' )

	for i in range(1,len(data)):
		truePos = 0.0
		trueNeg = 0.0
		falsePos = 0.0
		falseNeg = 0.0
		for j in range(i):
			if data[j].true == '1':
				truePos += 1.0
			if data[j].true == '-1':
				falsePos += 1.0
		for j in range(i,len(data)):
			if data[j].true == '-1':
				trueNeg += 1.0
			if data[j].true == '1':
				falseNeg += 1.0
		if trueNeg > 0 or falsePos > 0:
			rp1 = 1.0-(trueNeg/(trueNeg+falsePos))
		else:
			rp1 = 1.0
		if truePos > 0 or falseNeg > 0:
			rp2 = truePos/(truePos+falseNeg)
		else:
			rp2 = 0.0
		rocPoints.append( RocPoint(rp1, rp2 ) )

	for point in rocPoints:
		outFile.write(str(point.oneminusspec) + '\t' + str(point.sens) + '\n')
	
	sens = 0.0
	spec = 0.0
	for result in resultList:
		sens += result[0]
		spec += result[1]

	return (sens/len(resultList), spec/len(resultList), math.sqrt(spec/len(resultList)*sens/len(resultList)))


#################################################
# Run cross validation, build result files
#################################################

# Load data from positive and negative input files
allData = FeatureSet()
allData.load('data/'+posFile, patternClass='real')
allData.add_instances('data/'+negFile, patternClass='pseudo')
allData.libsvm_scale(paramOut = 'data/params')
subsets = allData.get_cv_subsets(numFolds)
resultList = []
# Go through all n folds...
for i in range(numFolds):
	# Build training and test sets
	testSet = subsets[i]
	trainSet = FeatureSet()
	for j in range(numFolds):
		if j != i:
			trainSet.add_instances_from_featureset(subsets[j])
	# Create svm files for train and test fold data. Train and test on these files.
	trainSet.weka_smote()
	trainSet.export_svm('data/trainSet.libsvm')
	testSet.export_svm('data/testSet.libsvm')
	# SVM settings for HMP features
	call('svm-train -c 1 -d 1 -h 1 -e 0.001 -g 0.06 -b 1 data/trainSet.libsvm models/'+str(i)+'.model', shell=True)
	# SVM settings for MicroPred features
	# call('svm-train -c 10000000 -d 1 -h 1 -e 0.001 -g 0.0019531 -b 1 data/trainSet.libsvm models/'+str(i)+'.model', shell=True)
	call('svm-predict -b 1 data/testSet.libsvm models/'+str(i)+'.model data/'+str(i)+'.results', shell=True)
	# Calculate sensitivity and specificity for fold model
	with open('data/'+str(i)+'.results', 'r') as resultFile:
		trueNeg = 0.0
		truePos = 0.0
		falseNeg = 0.0
		falsePos = 0.0
		for i in range(testSet.get_numpos()):
			if resultFile.readline()[0] == '1':
				truePos += 1.0
			else:
				falseNeg += 1.0
		for i in range(testSet.get_numneg()):
			if resultFile.readline()[0] == '1':
				falsePos += 1.0
			else:
				trueNeg += 1.0
		resultList.append( (truePos/(truePos+falseNeg),trueNeg/(trueNeg+falsePos)) )


###################
# Report Results
###################
for i in range(len(resultList)):
	print "## SVM "+str(i)+" ##"
	print 'Sensitivity: '+str(resultList[i][0])
	print 'Specificity: '+str(resultList[i][1])
print 'average Sensitivity: '+str(sum([result[0] for result in resultList])/numFolds)
print 'average Specificity: '+str(sum([result[1] for result in resultList])/numFolds)
print 'Geometric mean: '+str(pow(sum([result[0] for result in resultList])/numFolds*sum([result[1] for result in resultList])/numFolds, 0.5))
print '\n## Holdout test:', test_against_holdout(hoSpec, numFolds, posFile), "##"

#############
# Clean up
#############
call('rm data/params', shell=True)
call('rm data/testSet.libsvm', shell=True)
call('rm data/trainSet.libsvm', shell=True)
for i in range(numFolds):
	call('rm models/'+str(i)+'.model', shell=True)
	call('rm data/'+str(i)+'.results', shell=True)