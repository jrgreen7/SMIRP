# used to concatenate microPred features with seven additional HuntMi features

import sys

input1 = sys.argv[1]
input2 = sys.argv[2]
output = sys.argv[3]

output = open(output, 'w')

f1 = open(input1, 'r')
f2 = open(input2, 'r')

file1, file2 = f1.readlines(), f2.readlines()

for i in range(len(file2)):
  if 1:
    write = file1[i].replace('\n', '') + ',' + file2[i].replace('  ', ',')
    output.write(write)
  else:
    print "Could not merge this line: ", file1[i]
