# Extract sequences from RNALfold result file and export to fasta format file
def testLine(line):
	if len(line.split()) < 3:
		print "Invalid RNAfold data:", line
		return False
	if not line.split()[2].isdigit():
		print "Invalid RNAfold data:", line
		return False
	for c in line.split()[0]:
		if c not in '().':
			print "Invalid RNAfold data:", line
			return False
	return True

def RNAL_to_fasta(inPath, outPath, onePer = False):
	inFile = open(inPath, 'r')
	out = open(outPath, 'w')
	for line in inFile:
		if line[0] in 'AGCUT':
			for loc in geneLocs:
				outString = geneName+'-'+loc[0]+'\n' + line[int(loc[0])-1:int(loc[0])-1+int(loc[1])]+'\n'
				out.write(outString)
		if line[0] in '.()':
				if testLine(line) == True:
					if not onePer:
						geneLocs.append((line.split()[2], str(len(line.split()[0]))))
					if onePer:
						if done == False:
							geneLocs.append((line.split()[2], str(len(line.split()[0]))))
							done = True
		if line[0] == '>':
			geneName = line.split()[0]
			geneLocs = []
			if onePer:
				done = False
	inFile.close()
	out.close()

# Extract matching sequences from blast result file and export to fasta format file
def blast_to_fasta(inPath, outPath):
	inFile = open('inPath', 'r')
	out = open('outPath', 'w')

	for line in inFile:
		if line.split() != []:
			# Pull the title of the query out
			if line.split()[0] == 'Query=':
				currQuery = line.split()[1]
				matchNum = 0
			# Get the query
			if line[0] == '>':
				currSubject = line[1:-1]
				matchNum += 1
				newMatch = True
			# 
			if line.split()[0] == 'Sbjct' and newMatch:
				match = line.split()[2].replace('-', '')
				matchStart = line.split()[1]
				newMatch = False
			# 
			elif line.split()[0] == 'Sbjct' and not newMatch:
				match += line.split()[2].replace('-', '')
				out.write('>' + currQuery + '-match-' + str(matchNum) + ' ' + currSubject + ' '  + matchStart + '\n')
				out.write(match + '\n')
	inFile.close()
	out.close()

def RNAL_to_RNA(inPath, outPath):
	inFile = open('inPath', 'r')
	out = open('outPath', 'w')
	inFile.close()
	out.close()

def arff_to_svm(arff_fp, svm_fp):
	outFile = open(svm_fp, 'w')
	inFile = open(arff_fp, 'r')
	for line in inFile:
		if line[0] in '@ \n':
			continue
		else:
			if line.split(',')[-1] == 'real\n':
				outFile.write('1.0')
			else:
				outFile.write('-1.0')
			itemNum = 1
			for item in line.split(',')[:-1]:
				outFile.write(' '+str(itemNum)+':'+item)
				itemNum += 1
			outFile.write('\n')
	inFile.close()
	outFile.close()

def features_to_svm(feat_fp, svm_fp, value=None):
	outFile = open(svm_fp, 'w')
	inFile = open(feat_fp, 'r')	
	for line in inFile:
		if value != None:
			outFile.write(value)
		i = 1
		for item in line.split():
			if i != 1 or value != None:
				outFile.write(' ')
			outFile.write(str(i) + ':' + item)
			i += 1
		outFile.write('\n')
	inFile.close()
	outFile.close()

def csv_to_svm(csv_fp, svm_fp):
	outFile = open(svm_fp, 'w')
	inFile = open(csv_fp, 'r')
	for line in inFile:
		if line[0] == '"':
			continue
		else:
			if line.split(',')[-1] == '"real"\n':
				outFile.write('1.0')
			else:
				outFile.write('-1.0')
			itemNum = 1
			for item in line.split(',')[:-1]:
				outFile.write(' '+str(itemNum)+':'+item)
				itemNum+=1
			outFile.write('\n')
	inFile.close()
	outFile.close()