from collections import namedtuple
import random

Sequence = namedtuple('Sequence', ['name', 'seq'])

class SequenceList():

	def __init__(self):
		self.seqList = []

	def load_fasta(self, inPath):
		with open(inPath, 'r') as inFile:
			for line in inFile:
				if line[0] == '>':
					title = line.strip()
				else:
					seq = line.strip()
					self.seqList.append(Sequence(title, seq))

	def remove_all_redundant(self):
		tempList = []
		for i in range(len(self.seqList)):
			redundant = False
			for j in range(len(self.seqList)):
				if i != j:
					if self.seqList[i].seq in self.seqList[j].seq:
						if len(self.seqList[i].seq) < len(self.seqList[j].seq):
							redundant = True
							break
						elif i > j:
							redundant = True
							break
			if redundant == False:
				tempList.append(self.seqList[i])
		print "--Removed", str( len(self.seqList)-len(tempList) ), "sequences from the list."
		self.seqList = tempList

	def remove_redundant(self):
		newLength = 0
		oldLength = 1
		while (newLength != oldLength):
			tempList = []
			if not (self.seqList[0].seq in self.seqList[1].seq and len(self.seqList[0].seq) < len(self.seqList[1].seq)):
				tempList.append(self.seqList[0])
			for i in range(1,len(self.seqList)-1):
				redundant = False
				otherSeq = self.seqList[i-1]
				if self.seqList[i].seq in otherSeq.seq and self.seqList[i].seq != otherSeq.seq:
						redundant = True
				otherSeq = self.seqList[i+1]
				if self.seqList[i].seq in otherSeq.seq:
					redundant = True
				if not redundant:
					tempList.append(self.seqList[i])
			if not (self.seqList[-1].seq in self.seqList[-2].seq and len(self.seqList[-1].seq) < len(self.seqList[-2].seq)):
				tempList.append(self.seqList[-1])
			newLength = len(tempList)
			oldLength = len(self.seqList)
			print "--Removed", str( len(self.seqList)-len(tempList) ), "sequences from the list."
			self.seqList = tempList

	def select_random(self, num):
		if num < len(self.seqList):	
			selected = random.sample(range(len(self.seqList)), num)
			tempList = [self.seqList[x] for x in selected]
			self.seqList = tempList

	def export_fasta(self, outPath):
		with open(outPath, 'w') as outFile:
			for seq in self.seqList:
				outFile.write(seq.name+'\n')
				outFile.write(seq.seq+'\n')