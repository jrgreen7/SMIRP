# Runs extraction of 'loops' feature: the cumulative size 
# of internal loops found in the secondary structure

import re, sys

input = sys.argv[1]
output = sys.argv[2]

f = open(input, 'r')
f2 = open(output, 'w')

file = f.readlines()

count_lines = -10

for line in file:
  count_lines += 1
  if re.search('>', line):
    count_lines = 0
  if count_lines == 2:
    cut1_flag = False
    cut2_flag = False
    structure = line.split(' ')[0]
    cut1, cut2 = 0, 0
    for i in range(len(structure)):
      if structure[i] != '.' and not cut1_flag:
        cut1 = i
        cut1_flag = True
    for i in range(len(structure)):
      if structure[-i-1] != '.' and not cut2:
        cut2 = -i-1
    if not cut2:
      cut2 = -2
    if cut2 == -1:
      structure = structure[cut1:]    
    else:
      structure = structure[cut1:cut2+1]
    length = len(structure)
    middle = length/2
    half_structure_1 = structure[:middle]
    half_structure_2 = structure[middle:]
    reversed_brackets_1 = half_structure_1.count(')')
    reversed_brackets_2 = half_structure_2.count('(')
    reversed_brackets_1 = float(reversed_brackets_1)/middle
    reversed_brackets_2 = float(reversed_brackets_2)/middle
    reversed_brackets = (reversed_brackets_1 + reversed_brackets_2)/2
    f2.write(str(reversed_brackets) + '\n')  



    
    

