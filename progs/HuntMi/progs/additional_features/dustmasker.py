# Runs extraction of 'dm' feature: a percentage of low complexity
# regions detected in the sequence using Dustmasker

import re, sys, os

input = sys.argv[1]
output = sys.argv[2]

print "Running Dustmasker..."

dustmasker_input = input
dustmasker_output = input + '.dustmasker.tmp'

cmd = "./dustmasker -in " + dustmasker_input + " -out " + dustmasker_output + " -outfmt fasta -level 25"
os.system(cmd)

f = open(dustmasker_output, 'r')
f2 = open(output, 'w')
f3 = open('tmp', 'w')
file = f.readlines()

count = 0
for line in file:
  if re.search('>', line):
    count += 1
    if count > 1:
      f3.write('\n')
    f3.write(line)
  else:
    line = line.replace('\n', '')
    f3.write(line)

f3.write('\n')
f3.close()
f4 = open('tmp', 'r')
file4 = f4.readlines()

for line in file4:
  if re.search('>', line):
    small, big = 0, 0
  else: 
    if len(line) > 2:
      for elem in line:
        if re.search('[a|t|g|c|n]', elem):
          small += 1
        if re.search('[A|T|G|C|N]', elem):
          big += 1
        if elem == '\n':
          if not big:
            big = 1
          out = str(float(small)/float(len(line))) + '\n'
          f2.write(out)


                
