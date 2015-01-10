from classes.ResultSet import *
from classes.FeatureSet import *
from classes.ScanningTestFramework import *
import classes.FileConversion as fc
import classes.FastaOperations as fo
import classes.FoldOperations
from classes.SequenceList import *
from classes.MirnaDB import *
import random
import numpy
from collections import namedtuple
import matplotlib
import matplotlib.pyplot as plt
import scipy
import scipy.stats
import datetime
print datetime.datetime.now().time()

stf = ScanningTestFramework()
# stf.getValidMirnaHairpins('data/hsa_rna_regions.fasta.nr.hairpins', 'data/hsa_mirna_hp.15.18.75.fasta')
# stf.paramScan("test_param_out", 87)
# stf.setParams(-15.0, 18, 75)
stf.MPCorrelationMirna()
# stf.crossValidate('hsa_mirna_hp.15.18.75.fasta.micropred', 'hsa_coding_30k.15.18.75.fasta.micropred', 10)
# 
# def plot_roc(xvals, yvals, title):
# 	plt.scatter(xvals, yvals)
# 	plt.xlabel('1-Specificity')
# 	plt.ylabel('Sensitivity')
# 	plt.axis([0, 1, 0, 1])
# 	plt.title(title)
# 	# plt.show()
# 	plt.savefig('plot_test.png')

# headlines = open("hsa_scanning_model_params_test_results_2.1.8_87").readlines()[3::3]
# row_labels = sorted(list(set([float(x.split()[2]) for x in headlines])), reverse=True)
# print row_labels
# column_labels = sorted(list(set([int(x.split()[6]) for x in headlines])))

# datalines = open("hsa_scanning_model_params_test_results_2.1.8_87").readlines()[4::3]
# data = []
# for j in range(len(datalines)/22):
# 	data.append([float(x.split()[1]) for x in datalines[j*22:j*22+22]])
# len_array = numpy.array(data)
# data = []
# for j in range(len(datalines)/22):
# 	data.append([float(x.split()[1]) for x in datalines[j*22:j*22+22]])
# offset_array = numpy.array(data)

# datalines2 = open("hsa_scanning_model_params_test_results_2.1.8_87").readlines()[5::3]
# data = []
# for i in range(25):
# 	data = []
# 	for line in datalines2[22*i:22+i*22]:
# 		data.append( (float(line.split()[1]), float(line.split()[2])) )
# 	data = sorted(data, key=lambda x: x[1])
# 	xvals = [d[0] for d in data]
# 	yvals = [d[1] for d in data]
# 	color = numpy.random.rand()
# 	plt.scatter(xvals, yvals, c=numpy.random.rand(3,1))

# # plt.scatter(xvals, yvals)
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.title("Precision / recall for varying hairpin params L=120")
# plt.axis([0.55, 0.95, 0, 0.0015])
# plt.show()

# fig, ax = plt.subplots()
# ax.pcolor(len_array, cmap=matplotlib.cm.jet_r)
# ax.set_xticks(numpy.arange(len_array.shape[1])+0.5, minor=False)
# ax.set_yticks(numpy.arange(len_array.shape[0])+0.5, minor=False)
# ax.set_xticklabels(row_labels, minor=False)
# ax.set_yticklabels(column_labels, minor=False)
# ax.set_ylabel("Minimum number of base pairs")
# ax.set_xlabel("Maximum allowed free energy")
# plt.title("Offset difference (miRNA and extracted hairpin) L=120")
# plt.show()

# plot_roc('a', 'Title')

##########################################
# Scan through MFE and BP parameters
##########################################

# stf = ScanningTestFramework('data/hsamir.fasta', 'data/hsa_rna_regions.fasta', 'data/hsa_mature.fasta')

##########################################################
# Test predictions from extracted positive hairpins
# to see how many miRNA are recovered
##########################################################

# results = open("hsa_hairpins_L100.fasta.predict").readlines()
# results = [x[:-1] for x in results]
# mirna = open("hsa_hairpins_L100.fasta").readlines()[0:-1:2]
# mirna = ["-".join(x.split("-")[0:-1]) for x in mirna]
# z = zip(mirna, results)

# d = {}
# for pair in z:
# 	d[pair[0]] = "-1"

# for pair in z:
# 	if pair[1] == '1':
# 		d[pair[0]] = '1'

# num = 0
# for key in d.keys():
# 	if d[key] == '1':
# 		num += 1
# print num
# print len(d.keys())
# print float(num)/float(len(d.keys()))
# print d

#######################################
# Build positive set based on
# miRNA which form hairpins
# using RNAfold
#######################################

# lines = open("data/extracted positive hairpins/hsa_hairpins_L100.fasta", 'r').readlines()
# lines = ["-".join(l.split('-')[0:-1]) for l in lines[0:-1:2]]
# for l in lines:
# 	print l

# print len(lines)
# lines = list(set(lines))
# print len(lines)

# mirlines = open("data/extracted positive hairpins/hsamir.fasta", 'r').readlines()
# mirout = open("data/extracted positive hairpins/hsamir_L100.fasta", 'w')

# for line in lines:
# 	for i in range(len(mirlines)):
# 		if line == mirlines[i].split(' ')[0]:
# 			mirout.write(mirlines[i])
# 			mirout.write(mirlines[i+1])

#######################################
# Build positive sets for 'cousins'
# arabidopsis experiment
#######################################

# sl = SequenceList()
# sl.load_fasta('data/smo_neg.fasta.nr.hairpins')
# sl.select_random(5000)
# sl.export_fasta('data/smo_negative.fasta')

# mirdb = MiRNADB()
# mirdb.load_from_file('mirbase_19.mirdb')
# mirdb.remove_species('Arabidopsis thaliana')
# mirdb.remove_species('Arabidopsis lyrata')
# mirdb.remove_species('Brassica napus')
# mirdb.remove_species('Brassica oleracea')
# mirdb.remove_species('Brassica rapa')
# mirdb.remove_species('Gossypium arboreum')
# mirdb.remove_species('Gossypium herbaceum')
# mirdb.remove_species('Gossypium hirsutum')
# mirdb.remove_species('Gossypium raimondii')
# mirdb.remove_species('Theobroma cacao')
# mirdb.remove_species('Carica papaya')
# mirdb.remove_species('Cucumis melo')
# mirdb.remove_species('Hevea brasiliensis')
# mirdb.remove_species('Manihot esculenta')
# mirdb.remove_species('Ricinus communis')
# mirdb.remove_species('Acacia auriculiformis')
# mirdb.remove_species('Arachis hypogaea')
# mirdb.remove_species('Acacia mangium')
# mirdb.remove_species('Glycine max')
# mirdb.remove_species('Glycine soja')
# mirdb.remove_species('Lotus japonicus')
# mirdb.remove_species('Medicago truncatula')
# mirdb.remove_species('Phaseolus vulgaris')
# mirdb.remove_species('Vigna unguiculata')
# mirdb.remove_species('Bruguiera cylindrica')
# mirdb.remove_species('Bruguiera gymnorhiza')
# mirdb.remove_species('Malus domestica')
# mirdb.remove_species('Citrus clementine')
# mirdb.remove_species('Citrus reticulata')
# mirdb.remove_species('Citrus sinensis')
# mirdb.remove_species('Citrus trifoliata')
# mirdb.remove_species('Populus euphratica')
# mirdb.remove_species('Populus trichocarpa')
# mirdb.remove_species('Cynara cardunculus')
# mirdb.remove_species('Helianthus annuus')
# mirdb.remove_species('Helianthus argophyllus')
# mirdb.remove_species('Helianthus ciliaris')
# mirdb.remove_species('Helianthus exilis')
# mirdb.remove_species('Helianthus paradoxus')
# mirdb.remove_species('Helianthus petiolaris')
# mirdb.remove_species('Helianthus tuberosus')
# mirdb.remove_species('Digitalis purpurea')
# mirdb.remove_species('Rehmannia glutinosa')
# mirdb.remove_species('Salvia sclarea')
# mirdb.remove_species('Aquilegia caerulea')
# mirdb.remove_species('Nicotiana tabacum')
# mirdb.remove_species('Solanum lycopersicum')
# mirdb.remove_species('Solanum tuberosum')
# mirdb.remove_species('Vitis vinifera')
# mirdb.remove_species('Aegilops tauschii')
# mirdb.remove_species('Brachypodium distachyon')
# mirdb.remove_species('Elaeis guineensis')
# mirdb.remove_species('Festuca arundinacea')
# mirdb.remove_species('Hordeum vulgare')
# mirdb.remove_species('Oryza sativa')
# mirdb.remove_species('Sorghum bicolor')
# mirdb.remove_species('Saccharum officinarum')
# mirdb.remove_species('Saccharum ssp.')
# mirdb.remove_species('Triticum aestivum')
# mirdb.remove_species('Triticum turgidum')
# mirdb.remove_species('Zea mays')
# mirdb.remove_species('Chlamydomonas reinhardtii')
# mirdb.remove_species('Picea abies')
# mirdb.remove_species('Pinus densata')
# mirdb.remove_species('Pinus taeda')
# mirdb.remove_species('Physcomitrella patens')
# mirdb.remove_species('Selaginella moellendorffii')
# mirdb.filter_by_cluster('mirbase_19.clus', inSpecies='Ectocarpus siliculosus')
# mirdb.filter_top_results(600)
# mirdb.export_fasta('data/ath_cousins_8.fasta')

# fout = open("data/dmel-all-nsplit.fasta", 'w')
# lines = open("data/tmp/dmel-all.fasta.rna").readlines()
# linenum = 1
# for i in [1,3,5]:
# 	nsplit = lines[i].split("NNNN")
# 	for chunk in nsplit:
# 		if len(chunk) >= 40:
# 			if len(chunk) <= 5000:
# 				fout.write(">"+str(linenum)+"\n")
# 				linenum += 1
# 				fout.write(chunk+"\n")
# 			else:
# 				for i in range(len(chunk) / 5000):
# 					fout.write(">"+str(linenum)+"\n")
# 					linenum += 1
# 					fout.write(chunk[5000*i:5000*i+5000]+"\n")
					



# ////////////////////////////////////////////
# Histogram of MFE on neg micropred, neg DME
# ////////////////////////////////////////////

# micropredLines = open('data/micropred_neg.rnafold', 'r').readlines()[2::3]
# dmeLines = open('data/dme_negative.rnafold', 'r').readlines()[2::3]

# mpMfes = [float(line.split()[-1].strip('()\n')) for line in micropredLines]
# dmeMfes = [float(line.split()[-1].strip('()\n')) for line in dmeLines]*40

# size = len(mpMfes)
# x = scipy.arange(size)
# y = mpMfes
# size2 = len(dmeMfes)
# x2 = scipy.arange(size2)
# y2 = dmeMfes


# h2 = plt.hist(y2, bins=20, color='b')
# h = plt.hist(y, bins=20, color='r')
# plt.show()

# ////////////////////////////////////////////
# Histogram of MFE1 on pos/neg
# ////////////////////////////////////////////

# dmeLines = open('../Downloads/hmp-dme.csv', 'r').readlines()
# mposLines = open('../Downloads/micropred-pos.csv', 'r').readlines()
# mnegLines = open('../Downloads/micropred-neg.csv', 'r').readlines()

# dmePos = []
# dmeNeg = []
# mpPos = []
# mpNeg = []

# for line in dmeLines:
# 	if line.split(',')[-1].strip() == 'miRNA':
# 		dmePos.append( float(line.split(',')[1]) )
# 	else:
# 		dmeNeg.append( float(line.split(',')[1]) )

# for line in mposLines[1:]:
# 	mpPos.append( float(line.split(',')[1]) )
# for line in mnegLines[1:]:
# 	mpNeg.append( float(line.split(',')[1]) )

# posVals = [float(line.split()[0]) for line in posLines[1::2]]
# negVals = [float(line.split()[0]) for line in negLines[1::2]]

# print float(sum(posVals))/len(posVals)
# print float(sum(negVals))/len(negVals)


# with open('mirdup_scanning.tsv', 'w') as outFile:
# 	for i in range(50):
# 		numPos = 0.0
# 		numNeg = 0.0
# 		threshold = float(i)/100
# 		for posVal in posVals:
# 			if posVal < threshold:
# 				numPos += 1.0
# 		for negVal in negVals:
# 			if negVal < threshold:
# 				numNeg += 1.0
# 		outFile.write(str(threshold)+'\t'+str(numPos/len(posVals))+'\t'+str(numNeg/len(negVals))+'\n')


# size = len(dmePos)
# x = scipy.arange(size)
# y = dmePos
# size2 = len(dmeNeg)
# x2 = scipy.arange(size2)
# y2 = dmeNeg

# h = plt.hist(y, bins=20, color='r')
# h2 = plt.hist(y2, bins=20, color='b')
# plt.xlim(-0.025,0)
# plt.show()

# dist = getattr(scipy.stats, 'gamma')
# param = dist.fit(y)
# print param
# pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])
# plt.plot(pdf_fitted, label='beta'+" mammal")
# param2 = dist.fit(y2)
# print param2
# pdf_fitted2 = dist.pdf(x2, *param2[:-2], loc=param2[-2], scale=param2[-1])
# plt.plot(pdf_fitted2, label='beta'+" flower")
# plt.xlim(-2,10)
# plt.legend(loc='upper right')
# plt.show()


# ////////////////////////////////////////////
# Histogram of mirdup on pos/neg
# ////////////////////////////////////////////

# posLines = open('data/aca_positive.fasta.mirdup', 'r').readlines()
# negLines = open('data/aca_negative.fasta.mirdup', 'r').readlines()

# posVals = [float(line.split()[0]) for line in posLines[1::2]]
# negVals = [float(line.split()[0]) for line in negLines[1::2]]

# print float(sum(posVals))/len(posVals)
# print float(sum(negVals))/len(negVals)


# with open('mirdup_scanning.tsv', 'w') as outFile:
# 	for i in range(50):
# 		numPos = 0.0
# 		numNeg = 0.0
# 		threshold = float(i)/100
# 		for posVal in posVals:
# 			if posVal < threshold:
# 				numPos += 1.0
# 		for negVal in negVals:
# 			if negVal < threshold:
# 				numNeg += 1.0
# 		outFile.write(str(threshold)+'\t'+str(numPos/len(posVals))+'\t'+str(numNeg/len(negVals))+'\n')


# size = len(negVals)
# x = scipy.arange(size)
# y = negVals
# h = plt.hist(y, bins=20, color='r')
# size2 = len(posVals)
# x2 = scipy.arange(size2)
# y2 = posVals
# h2 = plt.hist(y2, bins=20, color='b')
# plt.show()

# /////////////////////////////////////////////
# Determining AU content distribution of miRNA
# /////////////////////////////////////////////
# pseudoSeqs = []

# lines = open('snail0_final.fasta', 'r').readlines()
# headers = [line.strip() for line in lines[0::2]]
# seqs = [line.strip() for line in lines[1::2]]
# for i in range(len(headers)):
# 	if float(headers[i].split()[-2]) > 0.98:
# 		pseudoSeqs.append(seqs[i])
# 	else:
# 		break
# lines = open('snail1_final.fasta', 'r').readlines()
# headers = [line.strip() for line in lines[0::2]]
# seqs = [line.strip() for line in lines[1::2]]
# for i in range(len(headers)):
# 	if float(headers[i].split()[-2]) > 0.98:
# 		pseudoSeqs.append(seqs[i])
# 	else:
# 		break
# lines = open('snail2_final.fasta', 'r').readlines()
# headers = [line.strip() for line in lines[0::2]]
# seqs = [line.strip() for line in lines[1::2]]
# for i in range(len(headers)):
# 	if float(headers[i].split()[-2]) > 0.98:
# 		pseudoSeqs.append(seqs[i])
# 	else:
# 		break

# lines = open("mirbase_19.mirdb", 'r').readlines()
# headers = [line[1:4] for line in lines[0::3]]
# seqs = [line.strip() for line in lines[1::3]]

# numAUAUAU = 0
# for seq in pseudoSeqs:
# 	if 'AUAUAUAUAUAUAUAUAU' in seq:
# 		numAUAUAU += 1

# AUcontents = []
# pseudoAUcontents = []

# for seq in seqs:
# 	AUcontent = 0.0
# 	AUs = seq.count('A') + seq.count('U')
# 	AUcontent = float(AUs) / len(seq)
# 	AUcontents.append(AUcontent)

# print "Mean AU content in miRNA:", float(sum(AUcontents))/len(AUcontents)

# for seq in pseudoSeqs:
# 	AUcontent = 0.0
# 	AUs = seq.count('A') + seq.count('U')
# 	AUcontent = float(AUs) / len(seq)
# 	pseudoAUcontents.append(AUcontent)

# print "Mean AU content in snail:", float(sum(pseudoAUcontents))/len(pseudoAUcontents)

# threshold = 0.70
# numAbove = 0
# for AUcontent in AUcontents:
# 	if AUcontent >= threshold:
# 		numAbove += 1
# print "miRNA above threshold:",numAbove
# print "Total miRNA:",len(AUcontents)
# print "Percent miRNA above:",float(numAbove)/len(AUcontents)

# with open('mirna_AU_thresholds.tsv', 'w') as outFile:
# 	for i in range(60, 91):
# 		threshold = float(i)/100
# 		numAbove = 0
# 		for AUcontent in AUcontents:
# 			if AUcontent >= threshold:
# 				numAbove += 1
# 		outFile.write(str(threshold)+'\t'+str(numAbove)+'\t'+str(float(numAbove)/len(AUcontents))+'\n')

# with open('snail_AU_thresholds.tsv', 'w') as outFile:
# 	for i in range(60, 91):
# 		threshold = float(i)/100
# 		numAbove = 0
# 		for AUcontent in pseudoAUcontents:
# 			if AUcontent >= threshold:
# 				numAbove += 1
# 		outFile.write(str(threshold)+'\t'+str(numAbove)+'\t'+str(float(numAbove)/len(pseudoAUcontents))+'\n')

# numAbove = 0
# for AUcontent in pseudoAUcontents:
# 	if AUcontent >= threshold:
# 		numAbove += 1
# print "Snail above threshold:",numAbove
# print "Snail with AUAUAU:", numAUAUAU
# print "Total Snail:",len(pseudoAUcontents)
# print "Percent snail above:",float(numAbove)/len(pseudoAUcontents)

# size2 = len(pseudoAUcontents)
# x2 = scipy.arange(size2)
# y2 = pseudoAUcontents
# h2 = plt.hist(y2, bins=20, color='b')
# size = len(AUcontents)
# x = scipy.arange(size)
# y = AUcontents
# h = plt.hist(y, bins=20, color='r')
# plt.show()

# /////////////////////////////////////////////////////
# Determining average lengths of miRNA for each species
# /////////////////////////////////////////////////////

# lines = open("mirbase_19.mirdb", 'r').readlines()
# headers = [line[1:4] for line in lines[0::3]]
# seqs = [line.strip() for line in lines[1::3]]

# data = zip(headers, seqs)
	
# numMirna = {}
# totalMirnaLen = {}
# mirnaLens = {}

# for datum in data:
# 	if datum[0] not in numMirna.keys():
# 		numMirna[datum[0]] = 0
# 		totalMirnaLen[datum[0]] = 0
# 		mirnaLens[datum[0]] = []
# 	numMirna[datum[0]] += 1
# 	totalMirnaLen[datum[0]] += len(datum[1])
# 	mirnaLens[datum[0]].append(len(datum[1]))

# # mammals = ['hsa','cfa','eca','mdo','meu','sha','age','lla','sla','mmu','mne','pbi','ggo','ppa','ppy','ptr','ssy','lca','oan','cgr','rno','bta','oar','ssc']
# mammals = ['hsa']
# flowers = ['cca','han','har','hci','hex','hpa','hpe','htu','aly','ata','bna','bol','bra','cpa','cme','hbr','mes','rco','aau','ahy','amg','gma','gso','lja','mtr','pvu','vun']
# virii = ['bhv','bkv','blv','bpc','ebv','hbv','hcm','hhv','hiv','hsv','hvs','hvt','ilt','jcv','ksh','mcm','mcv','mdv','mgh','prv','rlc','rrv','sv4']

# mammalData = []
# for key in mammals:
# 	for datum in mirnaLens[key]:
# 		mammalData.append(datum)

# flowerData = []
# for key in flowers:
# 	for datum in mirnaLens[key]:
# 		flowerData.append(datum)

# viriiData = []
# for key in virii:
# 	for datum in mirnaLens[key]:
# 		viriiData.append(datum)

# data = []
# for key in numMirna.keys():
# 	data.append( (key, totalMirnaLen[key]/numMirna[key]) )

# data.sort(key = lambda datum: datum[1], reverse=True)

# se = SpeciesEncoder('mirna_species_code.txt')

# outFile = open('mirbase_lengths.txt', 'w')
# for line in open('mirbase_hierarchy.txt', 'r'):
# 	if line.strip()[0] == '>':
# 		outFile.write(line)
# 	else:
# 		specName = line.split('(')[0].strip()
# 		specCode = se.get_code(specName)
# 		depth = line.rfind('\t') + 1
# 		for i in range(depth):
# 			outFile.write('\t')
# 		outFile.write(str(totalMirnaLen[specCode]/numMirna[specCode])+'\n')

# data = []
# for line in open('animal_lengths.txt', 'r'):
# 	data.append(int(line.strip()))



# size2 = len(flowerData)
# x2 = scipy.arange(size2)
# y2 = flowerData
# # h2 = plt.hist(y2, bins=100, color='b')

# size3 = len(viriiData)
# x3 = scipy.arange(size3)
# y3 = viriiData
# # h3 = plt.hist(y3, bins=10, color='g')

# size = len(mammalData)
# x = scipy.arange(size)
# y = mammalData
# # h = plt.hist(y, bins=20, color='r')

# dist_names = ['beta']

# for dist_name in dist_names:
# 	dist = getattr(scipy.stats, dist_name)
# 	param = dist.fit(y)
# 	print param
# 	pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])*size
# 	plt.plot(pdf_fitted, label=dist_name+" mammal")
# 	param2 = dist.fit(y2)
# 	print param2
# 	pdf_fitted2 = dist.pdf(x2, *param2[:-2], loc=param2[-2], scale=param2[-1])*size
# 	plt.plot(pdf_fitted2, label=dist_name+" flower")
# 	param3 = dist.fit(y3)
# 	print param3
# 	pdf_fitted3 = dist.pdf(x3, *param3[:-2], loc=param3[-2], scale=param3[-1])*size
# 	plt.plot(pdf_fitted3, label=dist_name+" virii")
# 	plt.xlim(0,500)
# plt.legend(loc='upper right')
# plt.show()

# fs = FeatureSet()
# fs.load('data/selected.negativeset.fa.-21.features')
# fs.export('data/micropred_neg.libsvm')

# fo.convert_DNA_to_RNA('data/hsa.fasta', 'data/hsa_rna_regions.fasta')

# # Load up the mirbase database and clusters, export clustered features.
# mirdb = MiRNADB()
# mirdb.load_from_file('mirbase_19.mirdb')

# mirdb.remove_species('Anolis carolinensis')
# mirdb.filter_by_cluster('mirbase_19.clus', inSpecies='Xenopus tropicalis')
# mirdb.filter_random_results(691)
# mirdb.export_features('data/xtrclust_randompos')
# mirdb.export_fasta('data/xtrclust_randompos.fasta')
# mirdb.export_labels('data/'+outPath+'.labels')

# //////////////////////////////////////////////////////
# Script for running RNAfold L=70 to L=120 for a given fasta
# //////////////////////////////////////////////////////

# with open("ath_L_70_to_400_offsets.txt", 'w') as outO:
# 	with open("ath_L_70_to_400_sizes.txt", 'w') as outD:
# 		for i in range(70,400):
# 			call('python extract_hairpins.py -i ath_rna_regions.fasta -n 6 -l '+str(i), shell=True)

# 			with open("data/ath_rna_regions.fasta.nr.hairpins", 'r') as hairpinFile:
# 				with open("data/ath_rna_regions.fasta", 'r') as regionFile:
# 					hairpinLines = hairpinFile.readlines()
# 					regionLines = regionFile.readlines()
# 					offsets = []
# 					diffs = []
# 					for j in range(100):
# 						totalLenDifference = 0
# 						totalOffset = 0
# 						numLens = 0
# 						hairpinHeaders = hairpinLines[0::2]
# 						hairpinData = [x.strip() for x in hairpinLines[1::2]]
# 						hairpinZip = zip(hairpinHeaders, hairpinData)
# 						hairpinSample = random.sample(hairpinZip, len(hairpinZip)/2)
# 						regionHeaders = regionLines[0::2]
# 						regionData = [x.strip() for x in regionLines[1::2]]
# 						regionZip = zip(regionHeaders, regionData)
# 						for hLine in hairpinSample:
# 							for rLine in regionZip:
# 								if rLine[0].strip() in hLine[0].strip():
# 									# outFile.write(hLine[0])
# 									lenDifference = len(rLine[1])-len(hLine[1])-59
# 									totalLenDifference += abs(lenDifference)
# 									numLens += 1
# 									startOffset = rLine[1].find(hLine[1])-30
# 									endOffset = len(rLine[1])-len(hLine[1])-59-startOffset
# 									totalOffset += (abs(startOffset)+abs(endOffset))
# 									# outFile.write(str(startOffset)+':')
# 									# outFile.write(str(endOffset)+'\n')
# 						offsets.append(float(totalOffset)/float(numLens))
# 						diffs.append(float(totalLenDifference)/float(numLens))
# 					offsetArray = numpy.array(offsets)
# 					diffArray = numpy.array(diffs)
# 					print diffArray
# 					outD.write(str(i)+'\t'+ str(sum(diffs)/len(diffs))+'\t'+str(numpy.std(diffArray))+'\n')
# 					outO.write(str(i)+'\t'+ str(sum(offsets)/len(offsets))+'\t'+str(numpy.std(offsetArray))+'\n')

# //////////////////////////////////////////////////////
# import sys
# from Bio import Entrez, SeqIO
# from Bio.Blast import NCBIWWW
# Entrez.email = "rpeace@sce.carleton.ca"

# mirnaLines = open("data/celmir.fasta","r").readlines()
# regionLines = open("data/tmp/cel.fasta.rna", "r").readlines()

# for line in regionLines:
# 	if line[0] == '>':
# 		title = line
# 	if line[0] != '>':
# 		found = False
# 		for mirLine in mirnaLines:
# 			if mirLine[:-1] in line or line in mirLine[:-1]:
# 				found = True
# 				# print "Found a match for "+title
# 		if found == False:
# 			print "Could not find a match for "+title


#######################################################################

# mirdb = MiRNADB()
# mirdb.load_from_file('mirbase_19.mirdb')
# mirdb.filter_by_cluster('mirbase_19.clus', inSpecies='Taeniopygia guttata')
# mirdb.remove_species('Gallus gallus')
# mirdb.filter_top_results(691)
# mirdb.export_fasta('data/finch_mirna.fasta')
# mirdb.export_labels('data/'+outPath+'.labels')

#######################################################################

# lines = open('data/positive_animal.fasta', 'r').readlines()
# with open('data/huntmi_animal_pos.fasta', 'w') as outFile:
# 	numSeqs = len(lines) / 2
# 	outLines = random.sample(range(numSeqs), 691)
# 	for i in range(numSeqs):
# 		if i in outLines:
# 			outFile.write(lines[2*i])
# 			outFile.write(lines[2*i+1])


# fs = FeatureSet()
# fs.load('data/ath_positive.fasta.hmp20')
# fs.export('data/ath_positive.libsvm')
# fs = FeatureSet()
# fs.load('data/ath_negative.fasta.hmp20')
# fs.export('data/ath_negative.libsvm')

# inPath = 'ebv_negative.fasta'
# outPath20 = 'data/' + inPath + '.hmp20'

# with open('data/'+inPath+'.hmp', 'w') as newOut:
# 	with open('data/'+inPath+'.headers', 'w') as headerOut:
# 		with open('data/'+inPath+'.csv', 'r') as oldOut:
# 			for line in oldOut:
# 				if line[0] in '1234567890-':
# 					newOut.write(line)
# 				else:
# 					headerOut.write(line)

# featureIndices = [1,3,9,10,20,26,27,32,33,34,35,39,37,43,44,210,153,156,169,203]

# with open('data/'+inPath+'.hmp', 'r') as inFile:
# 	with open(outPath20, 'w') as outFile:
# 		for line in inFile:
# 			features = line.split(',')
# 			for index in featureIndices:
# 				outFile.write(features[index]+',')
# 			outFile.write(features[-1])

print datetime.datetime.now().time()