from collections import namedtuple


Result = namedtuple('Result', ['header', 'seq', 'result', 'conf'])

class ResultSet():

	def __init__(self):
		self.results = []

	def load(self, fastaPath, resultPath):
		with open(fastaPath, 'r') as inFasta:
			with open(resultPath, 'r') as inResults:
				for line in inResults:
					if line[0] in '01-':
						h = inFasta.readline()[:-1]
						s = inFasta.readline()[:-1]
						r = line.split()[0]
						c = line.split()[1]
						self.results.append(Result(h,s,r,c))

	def sort(self):
		self.results = sorted(self.results, key=lambda result: float(result.conf), reverse=True)

	def export(self, outPath):
		with open(outPath, 'w') as outFile:
			for result in self.results:
				outFile.write(result.header)
				outFile.write(' [ Prediction: '+result.result+', Confidence: '+result.conf+' ]\n')
				outFile.write(result.seq+'\n')