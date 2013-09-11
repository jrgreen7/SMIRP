# Runs extraction of 'orf' feature: the maximal length of the amino acid 
# string without stop codons found in the sequence in three reading frames

import re, sys

input = sys.argv[1]
output = sys.argv[2]

f = open(input, 'r')
f2 = open(output, 'w')

file = f.readlines()

codons =   {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTA':'L', 
            'CTT':'L', 'CTG':'L', 'CTC':'L', 'ATT':'I', 'ATC':'I', 
            'ATA':'I', 'ATG':'M', 'GTA':'V', 'GTT':'V', 'GTG':'V', 
            'GTC':'V', 'TCA':'S', 'TCT':'S', 'TCG':'S', 'TCC':'S', 
            'AGT':'S', 'AGC':'S', 'CCA':'P', 'CCT':'P', 'CCG':'P', 
            'CCC':'P', 'ACA':'T', 'ACT':'T', 'ACG':'T', 'ACC':'T', 
            'GCA':'A', 'GCT':'A', 'GCG':'A', 'GCC':'A', 'TAT':'Y', 
            'TAC':'Y', 'TAA':'*', 'TAG':'*', 'TGA':'*', 'CAT':'H', 
            'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 
            'GAG':'E', 'TGT':'C', 'TGC':'C', 'TGG':'W', 'CGA':'R', 
            'CGT':'R', 'CGG':'R', 'CGC':'R', 'AGA':'R', 'AGG':'R', 
            'GGA':'G', 'GGT':'G', 'GGG':'G', 'GGC':'G', 'xxx':'X'}


maximal = []

for line in file:
  if not re.search('>', line) and len(line) > 2:
    seq = line.replace('\n', '')
    seq = seq.upper()
    seq = seq.replace('U', 'T')
    j = 0
    list = []
    for frame in range(0, 3):
      list.append(j)
      j = 0
      for i in range(len(seq)):
        if i%3 == frame:
          try:
            codon = seq[i:(i+3)]
            j += 1
            if codons[codon] == '*': 
              list.append(j)
              j = 0
          except:
            list.append(j)
            j = 0 
    maximal.append(max(list)/float(len(seq)))
    write = str(max(list)/float(len(seq))) + '\n'
    f2.write(write)
