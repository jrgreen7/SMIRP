### Runs extraction of seven new HuntMi miRNA features

import sys, os

fasta_file = sys.argv[1]

rnafold_output = '../results/' + fasta_file + '/' + fasta_file + '.fold'
loops_output = '../results/' + fasta_file + '/' + fasta_file + '.loops'
dustmasker_output = '../results/' + fasta_file + '/' + fasta_file + '.dustmasker'
micropred_output = '../results/' + fasta_file + '/selected_21_micropred_features'
triplet_output = '../results/' + fasta_file + '/' + fasta_file + '.fold.triplet'
translate_output = '../results/' + fasta_file + '/' + fasta_file + '.translate'
all_features = '../results/' + fasta_file + '/' + fasta_file + '.all_features'
fasta = fasta_file
fasta_file = '../data/' + fasta_file


print "Folding sequences..."
cmd = "./RNAfold/RNAfold < " + fasta_file + " > " + rnafold_output
os.system(cmd)

cmd = "rm *.ps"
os.system(cmd)

print "Triplets data..."
cmd = "python additional_features/triplet.py " + rnafold_output + ' ' + triplet_output
os.system(cmd)

print "Loops data..."
cmd = "python additional_features/loops.py " + rnafold_output + ' ' + loops_output
os.system(cmd)

print "Dustmasker data..."
cmd = "python additional_features/dustmasker.py " + fasta_file  + ' ' + dustmasker_output
os.system(cmd)

print "Translation data..."
cmd = "python additional_features/translate.py " + fasta_file + ' ' + translate_output
os.system(cmd)

print "merging with selected 21 microPred features..."
cmd = "python additional_features/merge.py " + micropred_output + ' ' + triplet_output + ' tmp1'
os.system(cmd)
cmd = "python additional_features/merge.py tmp1 " + loops_output + ' tmp2'
os.system(cmd)
cmd = "python additional_features/merge.py tmp2 " + dustmasker_output + ' tmp3'
os.system(cmd)
cmd = "python additional_features/merge.py tmp3 " + translate_output + ' ' + all_features
os.system(cmd)

