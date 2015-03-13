from collections import namedtuple
import random

MiRNA = namedtuple('MiRNA', ['name', 'sequence', 'features'])

class PhyloTree:
	def __init__(self, inName):
		self.species = []
		self.name = inName
		self.children = []

class SpeciesEncoder:
	def __init__(self, inPath):
		self.code_to_name = {}
		self.name_to_code = {}
		for line in open(inPath, 'r'):
			items = line.split(':')
			self.name_to_code[items[1].strip()] = items[0].strip()
			self.code_to_name[items[0].strip()] = items[1].strip()

	def get_name(self, inCode):
		return self.code_to_name[inCode]

	def get_code(self, inName):
		return self.name_to_code[inName]

	def encode_list(self, nameList):
		encList = []
		for l in nameList:
			encList.append([])
			for name in l:
				encList[-1].append(self.get_code(name))
		return encList

class SeqClusters:
	def __init__(self, inPath):
		self.clusters = []
		self.lens = []
		for line in open(inPath, 'r'):
			if line[0] == '>':
				self.clusters.append([])
				self.lens.append([])
			else:
				self.clusters[-1].append(line.split()[2].strip('.'))
				self.lens[-1].append(int(line.split()[1][0:-3]))

	def get_mean_rep(self, cluster, lens):
		averageLen = sum(lens)/len(lens)
		diffs = []
		for length in lens:
			diffs.append(abs(length-averageLen))
		topIdx = 0
		bestDiff = diffs[0]
		for j in range(1,len(diffs)):
			if diffs[j] < bestDiff:
				bestDiff = diffs[j]
				topIdx = j
		return cluster[topIdx]

	def get_representatives(self, species = "", minClusSize = 1):
		if species != "":
			reps = []
			for i in range(len(self.clusters)):
				if len(self.clusters[i]) >= minClusSize:
					filteredCluster = []
					filteredLens = []
					for classNum in range(len(species)):
						for j in range(len(self.clusters[i])):
							if self.clusters[i][j].split('-')[0][1:] in species[classNum]:
								filteredCluster.append(self.clusters[i][j])
								filteredLens.append(self.lens[i][j])
						if len(filteredCluster) > 0:
							break
					reps.append(self.get_mean_rep(filteredCluster, filteredLens))
		else:
			reps = []
			print minClusSize
			for i in range(len(self.clusters)):
				print self.clusters[i]
				if len(self.clusters[i]) >= minClusSize:
					reps.append(self.get_mean_rep(self.clusters[i], self.lens[i]))
		return reps

class SpeciesHierarchy:
	root = PhyloTree("Root")

	def __init__(self):
		self.build_hierarchy("mirbase_hierarchy_21.txt")

	def print_tree(self):
		self.print_tree_recurse(self.root, 0)

	def print_tree_recurse(self, node, depth):
		print '\t'*depth + node.name
		for s in node.species:
			print '\t'*(depth+1) + s
		for c in node.children:
			print_tree(c, depth+1)

	def add_taxonomy_to_list(self, inNode, inList):
		for spec in inNode.species:
			inList.append(spec)
		for child in inNode.children:
			self.add_taxonomy_to_list(child, inList)

	def species_by_similarity(self, inName):
		lengths = [[inName],[],[],[],[],[],[],[],[],[]]
		self.species_by_similarity_recurse(self.root, inName, lengths)
		for i in range(len(lengths)):
			lengths[i] = list(set(lengths[i]))
		return lengths

	def species_by_similarity_recurse(self, node, inName, lengths):
		stuffBelow = []
		returnList = []
		returnVal = False
		if inName in node.species:
			self.add_taxonomy_to_list(node, lengths[1])
			lengths[1].remove(inName)
			return 2
		for child in node.children:
			if self.species_by_similarity_recurse(child, inName, lengths) > 0:
				num = self.species_by_similarity_recurse(child, inName, lengths)
				self.add_taxonomy_to_list(node, lengths[num])
				for i in range(num):
					for spec in lengths[i]:
						if spec in lengths[num]:
							lengths[num].remove(spec)
				return num+1
		return 0

	def build_hierarchy(self, inPath):

		currNode = self.root
		inFile = open(inPath, 'r')

		for line in inFile:
			if line[0] == '>':
				self.root.children.append(PhyloTree(line.strip()))
			if line[1] == '>':
				self.root.children[-1].children.append(PhyloTree(line.strip()))
			if line[2] == '>':
				self.root.children[-1].children[-1].children.append(PhyloTree(line.strip()))
			if line[3] == '>':
				self.root.children[-1].children[-1].children[-1].children.append(PhyloTree(line.strip()))
			if line[4] == '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children.append(PhyloTree(line.strip()))
			if line[5] == '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children.append(PhyloTree(line.strip()))
			if line[6] == '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children.append(PhyloTree(line.strip()))
			if line[7] == '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children.append(PhyloTree(line.strip()))
			if line[8] == '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children.append(PhyloTree(line.strip()))

			if line.count('\t') == 1 and line[1] != '>':
				self.root.children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 2 and line[2] != '>':
				self.root.children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 3 and line[3] != '>':
				self.root.children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 4 and line[4] != '>':
				self.root.children[-1].children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 5 and line[5] != '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 6 and line[6] != '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 7 and line[7] != '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 8 and line[8] != '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())
			if line.count('\t') == 9 and line[9] != '>':
				self.root.children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].children[-1].species.append(line.split('(')[0].strip())

def filter_fasta(svmPath, labelsPath, repsPath, outPath):
	out = open(outPath, 'w')
	svmLines = open(svmPath, 'r').readlines()
	labels = [line.strip() for line in open(labelsPath, 'r').readlines()]
	reps = [line.strip() for line in open(repsPath, 'r').readlines()]
	for i in range(len(reps)):
		for j in range(len(labels)):
			if labels[j] == reps[i]:
				out.write(svmLines[j])
	out.close()

class MiRNADB:
	mirna_list = []
	encoder = SpeciesEncoder('mirna_species_code_21.txt')
	hierarchy = SpeciesHierarchy()

	def filter_random_results(self, numResults):
		self.mirna_list = random.sample(self.mirna_list, numResults)

	def filter_top_results(self, numResults):
		self.mirna_list = self.mirna_list[:numResults]

	def filter_to_species(self, inSpecies, inNum):
		specList = self.encoder.encode_list(self.hierarchy.species_by_similarity(inSpecies))
		new_mirna_list = []
		for specSubList in specList:
			print specSubList
			for mir in self.mirna_list:
				if mir.name[1:4] in specSubList:
					new_mirna_list.append(mir)
		new_mirna_list = new_mirna_list[:inNum]
		self.mirna_list = new_mirna_list

	def add_mirna(self, inName, inSequence, inFeatures):
		for mir in self.mirna_list:
			if mir.name == inName:
				return False
		self.mirna_list.append( MiRNA(inName, inSequence, inFeatures) )
		return True

	def remove_species(self, inSpecies):
		origLen = len(self.mirna_list)
		specCode = self.encoder.get_code(inSpecies)
		self.mirna_list = [mir for mir in self.mirna_list if not specCode in mir.name]
		print "Removed", origLen - len(self.mirna_list), "items from the database."

	# def filter_to_species(self, inSpecies):
	# 	specCode = self.encoder.get_code(inSpecies)
	# 	self.mirna_list = [mir for mir in self.mirna_list if specCode in mir.name]

	def filter_by_cluster(self, inClusPath, minClusSize=1, inSpecies=""):
		c = SeqClusters(inClusPath)
		if inSpecies == "":
			reps = c.get_representatives(minClusSize=minClusSize)
		else:
			l = self.hierarchy.species_by_similarity(inSpecies)
			l = self.encoder.encode_list(l)
			reps = c.get_representatives(minClusSize=minClusSize, species=l)
		tempList = []
		for rep in reps:
			for mir in self.mirna_list:
				if mir.name == rep:
					tempList.append(mir)
		self.mirna_list = []
		for mir in tempList:
			self.mirna_list.append(mir)
		self.mirna_list = [mir for mir in self.mirna_list if mir.name in reps]

	def load_from_file(self, inPath):
		with open(inPath, 'r') as inFile:
			name = "not None"
			while True:
				name = inFile.readline()
				if not name:
					break
				name = name.strip()
				seq = inFile.readline().strip()
				features = inFile.readline().strip()
				self.add_mirna(name, seq, features)

	def save_to_file(self, outPath):
		with open(outPath, 'w') as outFile:
			for mir in self.mirna_list:
				outFile.write(mir.name+'\n')
				outFile.write(mir.sequence+'\n')
				outFile.write(mir.features+'\n')

	def export_fasta(self, outPath):
		with open(outPath, 'w') as outFile:
			for mir in self.mirna_list:
				outFile.write(mir.name+'\n')
				outFile.write(mir.sequence+'\n')

	def export_labels(self, outPath):
		with open(outPath, 'w') as outFile:
			for mir in self.mirna_list:
				outFile.write(mir.name+'\n')

	def export_features(self, outPath):
		with open(outPath, 'w') as outFile:
			for mir in self.mirna_list:
				outFile.write(mir.features+'\n')

# print "Starting..."
# sh = SpeciesHierarchy()
# se = SpeciesEncoder('mirna_species_code.txt')
# sbs = sh.species_by_similarity('Homo sapiens')
# finalList = []
# for x in sbs:
# 	for y in x:
# 		finalList.append(se.get_code(y))
# print finalList

# # Opens a database, filters out stuff
# db = MiRNADB()
# db.load_from_file('mirbase_19.mirdb')
# db.remove_species('Rhesus lymphocryptovirus')
# db.filter_by_cluster('mirbase_19.clus', 3, 'Epstein Barr virus')
# db.filter_to_species('Epstein Barr virus', 691)
# db.export_fasta('mirbase_19_no_rlc_ebvclust.fasta')
# db.export_features('selected.mirbase_19_no_rlc_ebvclust.fasta.-21.features')
# print "Done!"
