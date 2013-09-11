import string


def filter_hairpins(inPath, outPath):
	with open(inPath, 'r') as fileIn:
		with open(outPath, 'w') as out:
			writeSeq = False
			for line in fileIn:
				if line[0] == '>':
					currHeader = line
					out.write(currHeader)
				elif line[0] in '.()':
					fold = line.split()[0]
					try:
					    float(line.split('(')[-1].split(')')[0])
					except ValueError:
					    print "Invalid float:", line.split('(')[-1].split(')')[0]
					    continue
					if '..' in line.split('(')[-1].split(')')[0]:
						print line
					# Restrictive Must-be-a-perfect-hairpin
					elif float(line.split('(')[-1].split(')')[0]) <= -15.00 and line.count('(') > 18:
						writeFold = True
						for s in fold.split(')')[1:]:
							if '(' in s:
								writeFold = False
						if writeFold:
							out.write(line)
							writeSeq = True
					# Less restrictive can-have-anything-inside-hairpin
					# elif float(line.split('(')[-1].split(')')[0]) <= -15.00  and fold.split(')')[0].count('(') >= 18 and fold.split('(')[-1].count(')') >= 18:
					# 	out.write(line)
					# 	writeSeq = True
				elif line[0] in 'AGCU' and writeSeq:
					#out.write(currHeader)
					out.write(line)
					writeSeq = False