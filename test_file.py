from subprocess import call
import sys, getopt

# Parameters:
# -i: Input fasta file
#
# Output is 
opts, extraparams = getopt.getopt(sys.argv[1:], 'i:')
for o,p in opts:
	if o == '-i':
		inPath = p

f = open(inPath, 'r')

lineNum = 0
for line in f:
	lineNum += 1
	if line[0] in 'ACGTUN':
		for c in line[:-1]:
			if c not in 'ACTGUN':
				print "Invalid character on line", str(lineNum), ":", c
	elif line[0] in '.()':
		for c in line[:-1]:
			if c not in '.()- 1234567890':
				print "Invalid character on line", str(lineNum), ":", c