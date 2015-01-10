from collections import namedtuple

Datum = namedtuple('Datum', ['name', 'conf', 'pClass'])

class PredictionData:

	def __init__(self):
		self.data = []
		self.rocData = []
		self.gmData = []
		return

	def load(self, inPath, pClass='1'):
		if inPath.split('.')[-1] == 'mirdup':
			self.load_mirdup(inPath, pClass)
		return

	def sort_by_confidence(self):
		self.data.sort(key=lambda datum:datum.conf, reverse=True)

	def build_roc_data(self):
		sortedData = sorted(self.data, key=lambda datum:datum.conf, reverse=True)
		self.rocData = []
		for i in range(len(sortedData)):
			TP = 0.0
			TN = 0.0
			FP = 0.0
			FN = 0.0
			for j in range(len(sortedData)):
				if j <= i and sortedData[j].pClass == '1':
					TP += 1
				if j > i and sortedData[j].pClass == '1':
					FN += 1
				if j <= i and sortedData[j].pClass == '-1':
					FP += 1
				if j > i and sortedData[j].pClass == '-1':
					TN += 1
			sens = TP / (TP + FN)
			spec = TN / (TN + FP)
			self.rocData.append((1-spec,sens))
			self.gmData.append(spec*sens)

	def load_mirdup(self, inPath, pClass='1'):
		with open(inPath, 'r') as inFile:
			for line in inFile:
				if line[0] == '>':
					datumName = line
				else:
					datumConf = float(line.strip().split()[0])
					datumClass = pClass
					self.data.append(Datum(datumName, datumConf, datumClass))

	def print_data(self):
		for datum in self.data:
			print datum.name,
			print str(datum.conf) + "\t" + datum.pClass

	def print_roc_data(self):
		self.build_roc_data()
		for datum in self.rocData:
			print str(datum[0]) + '\t' + str(datum[1])

	def print_gm_data(self):
		self.build_roc_data()
		for datum in self.gmData:
			print datum
		print "Max GM:", max(self.gmData)

p = PredictionData()
p.load('../data/aca_negative.fasta.mirdup', '-1')
p.load('../data/aca_positive.fasta.mirdup', '1')
p.sort_by_confidence()
p.print_data()