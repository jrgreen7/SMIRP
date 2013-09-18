import re, sys

fasta = sys.argv[1]

micropred = '../results/' + fasta + '/selected_21_micropred_features'

output_filename = '../results/' + fasta + '.removed'

arff = '../results/' + fasta + '.arff'

arff_filtered = '../results/' + fasta + '.filtered.arff'

fasta = '../data/' + fasta


f = open(micropred, 'r')
f2 = open(fasta, 'r')
out = open(output_filename, 'w')
filtered = open(arff_filtered, 'w')
arff = open(arff, 'r')

file = f.readlines()
file2 = f2.readlines()
file_arff = arff.readlines()

flag = False

count_lines = -1
remove = []

for line in file:
  count_lines += 1
  if re.search('\?,', line):
    remove.append(count_lines)

flag = False

count_lines = -1

for line in file2:
  if re.search('>', line):
    count_lines += 1
  if count_lines in remove:
    flag = True
  if count_lines not in remove:
    flag = False
  if flag:
    out.write(line)


for line in file_arff:
  if re.search('\?,', line):
    pass
  else:
    filtered.write(line)
    

