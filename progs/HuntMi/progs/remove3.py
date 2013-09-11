import re, sys

arff = sys.argv[1]

arff_filtered = arff.replace('.arff', '.filtered.arff')

filtered = open(arff_filtered, 'w')
arff = open(arff, 'r')

file_arff = arff.readlines()


for line in file_arff:
  if re.search('\?,', line):
    pass
  else:
    filtered.write(line)
    


