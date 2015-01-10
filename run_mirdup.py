from subprocess import call
import threading
import sys
import getopt
import os.path
import os

# Parameters:
# -i: Input fasta file
opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		inPath = p

inPath = 'data/snail_high_conf.fasta'
outPath = inPath+'.mirdup'

with open(outPath, 'w') as outFile:
	for line in open(inPath, 'r'):
		if line[0] == '>':
			outFile.write(line)
		else:
			seq = line.strip()
			os.chdir('progs/mirdup')
			call('java -Xms1500m -Xmx2500m -jar miRdup.jar -p -predict -u '+seq+' -d all.model -f out_tmp -r ViennaRNA-2.1.6/Progs/', shell=True)
			resultLines = open('out_tmp.generatedmirnas.all.model.miRdup.aln.txt', 'r').readlines()
			pos = float(resultLines[2].split()[3])
			total = float(resultLines[1].split()[3])
			best = float(resultLines[3].split()[1])
			outFile.write(str(pos/total)+'\t'+str(best)+'\n')
			os.chdir('../../')
