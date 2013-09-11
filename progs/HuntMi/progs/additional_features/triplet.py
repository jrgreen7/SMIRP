### Calculates the triplet frequencies in hairpin secondary structure.
### The triplets are: (((A, (((U, (((G, (((C
### (((A and )))A are the same

import re, sys

input = sys.argv[1]
output = sys.argv[2]

fold_file = input

f = open(fold_file, 'r')
f2 = open(output, 'w')

file = f.readlines()

dict_triplets = {}

i = -5
nr_of_elems = 0

for line in file:
  i += 1
  if re.search('>', line):
    i = 0
    dict_triplets = {}
    nr_of_elems = 0
  if i == 1:
    seq = line.replace('\n', '')
  if i == 2:
    structure = line.split(' ')[0]
    for j in range(len(structure) - 1):
      triplet = structure[j-1:j+2]
      triplet = triplet.replace(')', '(')
      nt = seq[j]
      elem = triplet + nt
      nr_of_elems += 1
      try:
        dict_triplets[elem] += 1
      except:
        dict_triplets[elem] = 1
    try:
      A = dict_triplets['(((A']/float(nr_of_elems)
    except:
      A = 0
    try:
      U = dict_triplets['(((U']/float(nr_of_elems)
    except:
      U = 0
    try:
      G = dict_triplets['(((G']/float(nr_of_elems)
    except:
      G = 0
    try:
      C = dict_triplets['(((C']/float(nr_of_elems)
    except:
      C = 0
    write = str(A) + '  ' + str(U) + '  ' + str(G) + '  ' + str(C) + '\n'
    f2.write(write)


