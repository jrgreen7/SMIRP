import sys, os


os.chdir("progs")

def check_file(filename):
  ''' checks if file exists '''
  if(os.path.exists(filename)):
    ''' checks the file size '''
    filesize = os.path.getsize(filename)
  else:
    filesize = 0
  return filesize


if len(sys.argv) == 3:

  positive = sys.argv[1]
  negative = sys.argv[2]

  if check_file('../data/' + positive) and check_file('../data/' + negative):
    pass
  else:
    info = 'input files \"' + positive + '\" and \"' + negative + '\" not found at ../data/! Had to abort'
    sys.exit(info)


  arff_file = '../results/' + positive + '-' + negative + '.arff'

  positive_features = '../results/' + positive + '/' + positive + '.all_features'
  negative_features = '../results/' + negative + '/' + negative + '.all_features'


  print "\nCALCULATING microPred FEATURES FOR ", positive
  cmd = "python micropred_features.py " + positive 
  os.system(cmd)

  print "\nCALCULATING microPred FEATURES FOR ", negative
  cmd = "python micropred_features.py " + negative
  os.system(cmd)

  print "\nCALCULATING ADDITIONAL FEATURES FOR ", positive
  cmd = "python additional_features.py " + positive
  os.system(cmd)

  print "\nCALCULATING ADDITIONAL FEATURES FOR ", negative
  cmd = "python additional_features.py " + negative
  os.system(cmd)

  print "\nCREATING .arff FILE"
  cmd = "python additional_features/arff.py " + positive_features + ' ' + negative_features + ' ' + arff_file
  os.system(cmd)

  cmd = "rm " + positive + '.*'
  os.system(cmd)
  cmd = "rm " + negative + '.*'
  os.system(cmd)

  if check_file(arff_file):
    print "Results saved at ", arff_file
  else:
    print "Failed to generate output file!"

  os.system('python remove2.py ' + positive)
  os.system(cmd)

  os.system('python remove2.py ' + negative)
  os.system(cmd)

  os.system('python remove3.py ' + arff_file)
  os.system(cmd)







if len(sys.argv) == 2:

  input = sys.argv[1]

  if check_file('../data/' + input):
    pass
  else:
    info = 'input file \"' + input + '\" not found at ../data/! Had to abort'
    sys.exit(info)

  arff_file = '../results/' + input + '.arff'

  input_features = '../results/' + input + '/' + input + '.all_features'


  print "\nCALCULATING microPred FEATURES FOR ", input 
  cmd = "python micropred_features.py " + input
  os.system(cmd)

  print "\nCALCULATING ADDITIONAL FEATURES FOR ", input 
  cmd = "python additional_features.py " + input 
  os.system(cmd)

  print "\nCREATING .arff FILE"
  cmd = "python additional_features/arff.py " + input_features + ' ' + arff_file
  os.system(cmd)

  cmd = "rm " + input + '.*'
  os.system(cmd)

  if check_file(arff_file):
    print "Results saved at ", arff_file
  else:
    print "Failed to generate output file!"

  # os.system('python remove.py ' + input)
  # os.system(cmd)



if len(sys.argv) == 1:
  print "No input filenames provided!"
  print "USAGE: \npython extractFeatures.py filename"
  print "\tOR\npython extractFeatures.py positive_filename negative_filename"

if len(sys.argv) > 3:
  print "Too many arguments provided"
  print "USAGE: \npython extractFeatures.py filename"
  print "\tOR\npython extractFeatures.py positive_filename negative_filename"


os.system('rm ../data/*.tmp')
os.system('rm tmp*')


