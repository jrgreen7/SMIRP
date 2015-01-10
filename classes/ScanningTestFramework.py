from collections import namedtuple
from subprocess import call
import string
import random
import math
from collections import namedtuple
from operator import attrgetter
from classes.FeatureSet import FeatureSet
from classes.Plotter import Plotter

HairpinMetrics = namedtuple('HairpinMetrics', ['len', 'offset'])
MirnaSeqs = namedtuple('MirnaSeqs', ['five', 'three', 'pre', 'region'])
Result = namedtuple('Result', ['t', 'p', 'conf'])

class ScanningTestFramework:

	def __init__(self, mirnaPath="", regionPath="", maturePath=""):
		if mirnaPath == "":
			mirnaPath = 'data/hsamir.fasta'
			regionPath = 'data/hsa_rna_regions.fasta'
			maturePath = 'data/hsa_mature.fasta'

		self.mirna = {}
		mirnaZip = self.zipFasta(mirnaPath)
		regionZip = self.zipFasta(regionPath)
		matureZip = self.zipFasta(maturePath)

		self.ci = 0.0
		self.hpSens = 0.0
		self.hpSpec = 0.0

		for m in mirnaZip:
			self.mirna[m[0].split()[0]] = MirnaSeqs('','',m[1],'')

		for m in regionZip:
			self.mirna[m[0]] = self.mirna[m[0]]._replace(region = m[1])

		for mirName in self.mirna.keys():
			if len(mirName.split('-')) == 4:
				baseMirName = '-'.join(mirName.split('-')[:-1])
				for mat in matureZip:
					matName = mat[0].split()[0]
					if baseMirName == matName.strip('-5p').strip('-3p'):
						if '-5p' in matName:
							self.mirna[mirName] = self.mirna[mirName]._replace(five = mat[1])
						elif '-3p' in matName:
							self.mirna[mirName] = self.mirna[mirName]._replace(three = mat[1])
						else:
							self.mirna[mirName] = self.mirna[mirName]._replace(five = mat[1])
							self.mirna[mirName] = self.mirna[mirName]._replace(three = mat[1])
			for mat in matureZip:
				if mirName == matName.strip('-5p').strip('-3p'):
					if '-5p' in matName:
						self.mirna[mirName] = self.mirna[mirName]._replace(five = mat[1])
					elif '-3p' in matName:
						self.mirna[mirName] = self.mirna[mirName]._replace(three = mat[1])
					else:
						self.mirna[mirName] = self.mirna[mirName]._replace(five = mat[1])
						self.mirna[mirName] = self.mirna[mirName]._replace(three = mat[1])

		self.flankSize = 30

	def setParams(self, m, p, hpLen):
		call('python extract_hairpins.py -i hsa_coding_30000-35000.fasta -n 6 -l '+str(hpLen)+' -m 0.00 -p 1', shell=True)
		codingAll = self.zipFasta('data/hsa_coding_30000-35000.fasta.nr.hairpins')
		call('python extract_hairpins.py -i hsa_coding_30000-35000.fasta -n 6 -l '+str(hpLen)+' -m '+str(m)+' -p '+str(p), shell=True)
		codingSome = self.zipFasta('data/hsa_coding_30000-35000.fasta.nr.hairpins')
		call('python extract_hairpins.py -i hsa_rna_regions.fasta -n 6 -l '+str(hpLen)+' -m '+str(m)+' -p '+str(p), shell=True)
		self.getValidMirnaHairpins('data/hsa_rna_regions.fasta.nr.hairpins', 'data/hsa_rna_regions.fasta.nr.validhp')
		posZip = self.zipFasta('data/hsa_rna_regions.fasta.nr.validhp')
		self.hpSens = float(len(posZip)) / float(len(self.mirna.keys()))
		self.hpSpec = float(len(codingSome)) / float(len(codingAll))
		call('python extract_hairpins.py -i hsa_chr11_sample_2.fasta -n 6 -l '+str(hpLen)+' -m '+str(m)+' -p '+str(p), shell=True)
		negZip = self.zipFasta('data/hsa_chr11_sample_2.fasta.nr.hairpins')
		numPos = len(self.mirna.keys())*2
		numNeg = len(negZip) * 10000
		print 'Number of microRNA (estimated): '+str(numPos)
		print 'Number of hairpins in genome (estimated): '+str(numNeg)
		self.ci = float(numNeg)/float(numPos)
		print 'Class imbalance (estimated): '+str(self.ci)
		print 'Sensitivity of hairpin extraction :'+str(self.hpSens)
		print 'Specificity of hairpin extraction :'+str(self.hpSpec)

	def paramScan(self, outPath, hpLen):
		outFile = open(outPath, 'w')
		outFile.write('Base: (Min MFE: 0.00 Min paired bases: 1)\n')
		call('python extract_hairpins.py -i hsa_rna_regions.fasta -n 6 -l '+str(hpLen)+' -m 0.00 -p 1', shell=True)
		metrics = self.getHairpinMetrics('data/hsa_rna_regions.fasta.nr.hairpins')
		outFile.write(str(metrics.len)+"\t")
		outFile.write(str(metrics.offset)+"\n")
		outFile.write(str(self.getHairpinSensitivity('data/hsa_rna_regions.fasta.nr.hairpins'))+"\t")
		call('python extract_hairpins.py -i hsa_chr11_sample_2.fasta -n 6 -l '+str(hpLen)+' -m 0.00 -p 1', shell=True)
		negLines = open('data/hsa_chr11_sample_2.fasta.nr.hairpins', 'r').readlines()
		outFile.write(str(len(negLines)/2)+"\n")
		for i in range(5,30):
			for j in [-7.0, -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, -14.0, -15.0, -16.0, -17.0, -18.0, -19.0, -20.0, -21.0, -22.0, -23.0, -24.0, -25.0, -26.0, -27.0, -28.0]:
				outFile.write('Min MFE: '+str(j)+' Min paired bases: '+str(i)+'\n')
				call('python extract_hairpins.py -i hsa_rna_regions.fasta -n 6 -l '+str(hpLen)+' -m '+str(j)+' -p '+str(i), shell=True)
				metrics = self.getHairpinMetrics('data/hsa_rna_regions.fasta.nr.hairpins')
				outFile.write(str(metrics.len)+"\t")
				outFile.write(str(metrics.offset)+"\n")
				call('python extract_hairpins.py -i hsa_chr11_sample_2.fasta -n 6 -l '+str(hpLen)+' -m '+str(j)+' -p '+str(i), shell=True)
				negLines = open('data/hsa_chr11_sample_2.fasta.nr.hairpins', 'r').readlines()
				outFile.write(str(len(negLines)/2)+"\t")
				outFile.write(str(self.getHairpinSensitivity('data/hsa_rna_regions.fasta.nr.hairpins'))+"\t")
				outFile.write(str(float(self.getHairpinSensitivity('data/hsa_rna_regions.fasta.nr.hairpins')*5000)/float(len(negLines)*5000))+'\n')
		outFile.close()

	def getValidMirnaHairpins(self, hairpinPath, outPath):
		hairpinZip = self.zipFasta(hairpinPath)
		outFile = open(outPath, 'w')

		for mirName in self.mirna.keys():
			for hairpin in hairpinZip:
				if '-'.join(hairpin[0].split('-')[:-1]) == mirName:
					if self.mirna[mirName].five in hairpin[1] and self.mirna[mirName].three in hairpin[1]:
						outFile.write('>'+hairpin[0]+'\n')
						outFile.write(hairpin[1]+'\n')
						break

	def zipFasta(self, path):
		lines = open(path, 'r').readlines()
		headers = [x.strip()[1:] for x in lines[0::2]]
		data = [x.strip() for x in lines[1::2]]
		return zip(headers, data)

	def getHairpinMetrics(self, hairpinPath):
		hairpinZip = self.zipFasta(hairpinPath)
		offsets = []
		diffs = []
		totalLenDiff = 0
		totalOffset = 0
		numLens = 0
		for hLine in hairpinZip:
			mirName = '-'.join(hLine[0].split('-')[:-1])
			region = self.mirna[mirName].region
			totalLenDiff += abs(len(region) - len(hLine[1]) - (self.flankSize*2-1))
			numLens += 1
			startOffset = region.find(hLine[1])-self.flankSize
			endOffset = len(region)-len(hLine[1]) - (self.flankSize*2-1)-startOffset
			totalOffset += (abs(startOffset) + abs(endOffset))

		metrics = HairpinMetrics(float(totalLenDiff)/float(numLens),float(totalOffset)/float(numLens))
		return metrics

	def getHairpinSensitivity(self, m, p, hpLen):
		hairpinZip = self.zipFasta(hairpinPath)
		numRecovered = 0
		numTotal = len(self.mirna.keys())

		for mirName in self.mirna.keys():
			for hairpin in hairpinZip:
				if '-'.join(hairpin[0].split('-')[:-1]) == mirName:
					if self.mirna[mirName].five in hairpin[1] and self.mirna[mirName].three in hairpin[1]:
						numRecovered += 1
						break

		return float(numRecovered) / float(numTotal)

	def getHairpinSensitivity(self, hairpinPath):
		hairpinZip = self.zipFasta(hairpinPath)
		numRecovered = 0
		numTotal = len(self.mirna.keys())

		for mirName in self.mirna.keys():
			for hairpin in hairpinZip:
				if '-'.join(hairpin[0].split('-')[:-1]) == mirName:
					if self.mirna[mirName].five in hairpin[1] and self.mirna[mirName].three in hairpin[1]:
						numRecovered += 1
						break

		return float(numRecovered) / float(numTotal)

	def crossValidate(self, posFile, negFile, numFolds):
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
				with open("data/"+str(i)+".sresults", 'w') as resultOut:
					# resultLines = resultFile.readlines()
					# posLines = resultLines[1:testSet.get_numpos())].sorted( key=lambda l: float(l.split()[1]) )
					# negLines = resultLines[testSet.get_numpos():].sorted( key=lambda l: float(l.split()[1]) )
					trueNeg = 0.0
					truePos = 0.0
					falseNeg = 0.0
					falsePos = 0.0
					resultSet = []
					resultFile.readline()
					for j in range(testSet.get_numpos()):
						line = resultFile.readline()
						if line[0] == '1':
							resultSet.append(Result(t='1', p='1', conf=line.split()[1]))
							truePos += 1.0
						else:
							resultSet.append(Result(t='1', p='0', conf=line.split()[1]))
							falseNeg += 1.0
					for j in range(testSet.get_numneg()):
						line = resultFile.readline()
						if line[0] == '1':
							resultSet.append(Result(t='0', p='1', conf=line.split()[1]))
							falsePos += 1.0
						else:
							resultSet.append(Result(t='0', p='0', conf=line.split()[1]))
							trueNeg += 1.0
					resultSet = sorted(resultSet, key=lambda l: float(l.conf), reverse=True)
					for r in resultSet:
						resultOut.write(r.t + '\t' + r.p + '\t' + r.conf + '\n')

					resultList.append( (truePos/(truePos+falseNeg),trueNeg/(trueNeg+falsePos)) )

					with open("roc_"+str(i)+".tsv", 'w') as rocOut:
						with open("pr_"+str(i)+".tsv", 'w') as prOut:
							ssList = []
							prList = []
							sens = 0.0
							spec = 1.0
							for r in resultSet:
								if r.t == '1':
									sens += 1.0 / testSet.get_numpos()
								if r.t == '0':
									spec -= 1.0 / testSet.get_numneg()
								ssList.append((sens*self.hpSens, (1-spec)*self.hpSpec))
								prList.append((sens*self.hpSens/(sens*self.hpSens+(1-spec)*self.ci*self.hpSpec), sens*self.hpSens))
								rocOut.write(str(sens)+'\t'+str(1-spec)+'\n')
								prOut.write(str(sens/(sens+spec*self.ci))+'\t'+str(sens)+'\n')

					p = Plotter()
					p.plot_roc(ssList, "Test", "roc_"+str(i)+".png")
					p.plot_pr(prList, "Test", self.ci, "pr_"+str(i)+".png")

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

	def MPCorrelationMirna(self):
		call('progs/ViennaRNA-2.1.8/Progs/RNAfold -d2 --noLP < data/hsa_mirna_hp.15.18.100.fasta > data/hsa_rna_regions.folds', shell=True)
		MPpairs = []
		with open('data/hsa_rna_regions.folds', 'r') as fileIn:
			with open('results/mp_pairs_pos.tsv', 'w') as fileOut:
				writeSeq = False
				for line in fileIn:
					if line[0] in '.()':
						writeFold = True
						fold = line.split()[0]
						try:
						    float(line.split('(')[-1].split(')')[0])
						except ValueError:
						    print "Invalid float:", line.split('(')[-1].split(')')[0]
						    continue
						if '..' in line.split('(')[-1].split(')')[0]:
							print line
						# Restrictive Must-be-a-perfect-hairpin
						for s in fold.split(')')[1:]:
							if '(' in s:
								writeFold = False
						if writeFold:
							MPpairs.append( (float(line.split('(')[-1].split(')')[0]), line.count('(') - 1) )
							fileOut.write( str(float(line.split('(')[-1].split(')')[0])) +'\t'+ str(line.count('(')-1)+'\n')
						# Less restrictive can-have-anything-inside-hairpin
						# elif float(line.split('(')[-1].split(')')[0]) <= -15.00  and fold.split(')')[0].count('(') >= 18 and fold.split('(')[-1].count(')') >= 18:
						# 	out.write(line)
						# 	writeSeq = True
		p = Plotter()
		p.plot_scatter(MPpairs, "mp_scatter_pos.png")


	# ###########################################
	# Draws a scatter plot of MFE and BP values for poops
	# ###########################################
	def MPCorrelation(self):
		call('progs/ViennaRNA-2.1.8/Progs/RNALfold -d2 --noLP -L 100 < data/hsa_chr11_sample_2.fasta > data/hsa_chr11_sample_2.folds', shell=True)
		MPpairs = []
		with open('data/hsa_chr11_sample_2.folds', 'r') as fileIn:
			with open('results/mp_pairs_neg.tsv', 'w') as fileOut:
				writeSeq = False
				for line in fileIn:
					if line[0] in '.()':
						writeFold = True
						fold = line.split()[0]
						try:
						    float(line.split('(')[-1].split(')')[0])
						except ValueError:
						    print "Invalid float:", line.split('(')[-1].split(')')[0]
						    continue
						if '..' in line.split('(')[-1].split(')')[0]:
							print line
						# Restrictive Must-be-a-perfect-hairpin
						for s in fold.split(')')[1:]:
							if '(' in s:
								writeFold = False
						if writeFold:
							MPpairs.append( (float(line.split('(')[-1].split(')')[0]), line.count('(') - 1) )
							fileOut.write( str(float(line.split('(')[-1].split(')')[0])) +'\t'+ str(line.count('(')-1)+'\n')
						# Less restrictive can-have-anything-inside-hairpin
						# elif float(line.split('(')[-1].split(')')[0]) <= -15.00  and fold.split(')')[0].count('(') >= 18 and fold.split('(')[-1].count(')') >= 18:
						# 	out.write(line)
						# 	writeSeq = True
		p = Plotter()
		p.plot_scatter(MPpairs, "mp_scatter_neg.png")
