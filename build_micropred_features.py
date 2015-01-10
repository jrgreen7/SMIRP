import threading
import sys
import getopt
from subprocess import call
import os.path
import classes.FastaOperations as FastaOps
import os

opts, extraparams = getopt.getopt(sys.argv[1:], 'i:n:l:')
for o,p in opts:
	if o == '-i':
		inPath = p
	if o == '-n':
		numThreads = int(p)


class myThread(threading.Thread):
	def __init__(self, inPath):
		threading.Thread.__init__(self)
		self.inPath = inPath

	def run(self):
		with open("../data/"+self.inPath, "r") as f:
			fastaLines = f.readlines()

			if os.path.exists("../data/selected."+self.inPath+".-21.features"):
				prevLines = open("../data/selected."+self.inPath+".-21.features").readlines()
			else:
				prevLines = []

			for i in range(1,len(fastaLines),2):
				header = fastaLines[i-1]
				seq = fastaLines[i]
				if header not in prevLines:
					with open("../data/tmp_"+self.inPath, 'w') as tempOut:
						tempOut.write(header+seq)
					
					call("perl microPred.pl tmp_"+self.inPath, shell=True)
					
					with open("../data/selected.tmp_"+self.inPath+".-21.features", "r") as inFile:
						with open("../data/selected."+self.inPath+".-21.features", 'a') as finalOut:
							features = inFile.readline()
							if features != '':
								finalOut.write(header)
								finalOut.write(features)

###########################################
# Main script starts here
###########################################

# Move data file into microPred directory
call('cp data/'+inPath+' progs/microPred/data/'+inPath, shell=True)
exitFlag = 0

# Split input fasta into <numThreads> smaller fasta files
FastaOps.split_fasta('progs/microPred/data/'+inPath, numThreads)

#Get a path to the input file that doesn't have '.fasta' on the end
threadPath = inPath.rsplit('.',1)[0]
threadExt = inPath.rsplit('.',1)[1]
# Move to microPred directory and run microPred on all the smaller files (multi-threaded)
os.chdir('progs/microPred/progs')
threads = []
for i in range(numThreads):
	threads.append(myThread(threadPath+'.'+str(i)+'.'+threadExt))

for thread in threads:
	thread.start()

for thread in threads:
	thread.join()
os.chdir('../../../')

# 
with open('progs/microPred/data/selected.'+inPath+'.-21.features', 'w') as finalOut:
	for i in range(numThreads):
		with open('progs/microPred/data/selected.'+threadPath+'.'+str(i)+'.'+threadExt+'.-21.features', 'r') as finalIn:
			for line in finalIn:
				finalOut.write(line)
		call('rm progs/microPred/data/selected.'+threadPath+'.'+str(i)+'.'+threadExt+'.-21.features', shell=True)

with open('data/'+inPath+'.micropred', 'w') as newOut:
	with open('data/'+inPath+'.headers', 'w') as headerOut:
		with open('progs/microPred/data/selected.'+inPath+'.-21.features', 'r') as oldOut:
			for line in oldOut:
				if line[0] in '1234567890-':
					newOut.write(line)
				else:
					headerOut.write(line)

# call('rm progs/microPred/data/*', shell=True)