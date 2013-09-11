import sys, os, fileinput, re, time

inputfile = sys.argv[1]

count = 0 #counts analyzed sequences
inputfile_files = inputfile + '.files'
inputfile_path = '../data/' + inputfile

cmd = 'mkdir ../results/' + inputfile
os.system(cmd)

cmd = 'touch ../results/' + inputfile + '/selected_21_micropred_features'
os.system(cmd)

output = '../results/' + inputfile + '/selected_21_micropred_features'
out_main = open(output, 'w')


#----------------------------------------------------------------------------------

for line in fileinput.input(inputfile_path):
  if re.search('>', line):
    count += 1
    print "Analysing sequence nr", count
  else:
    cmd = "rm ../data/" + inputfile_files + "*"
    os.system(cmd)
    outfile = '../data/' + inputfile_files
    out = open(outfile, 'w')
    write = '>seq\n' + line
    out.write(write)
    out.close()

    cmd = "./RNAfold/RNAfold < " + '../data/' + inputfile_files + " > " + '../data/' + inputfile_files + ".fold" 
    os.system(cmd)
    
    #----------------------------------------------------------------------------------
    
    #calculating miPred sequential and structural features
    cmd = "perl miPred/genRNAStats.pl < ../data/" + inputfile_files + " > ../data/" + inputfile_files + ".data1"
    os.system(cmd)
    
    cmd = "miPred/./RNAspectral.exe < ../data/" + inputfile_files + ".fold > ../data/" + inputfile_files + ".data2"
    os.system(cmd)
    
    #----------------------------------------------------------------------------------

    #calculating z-features
    cmd = "perl miPred/genRandomRNA.pl -n 1000 -m d < ../data/" + inputfile_files + " > ../data/" + inputfile_files + ".random.fasta"
    os.system(cmd)

    cmd = "./RNAfold/RNAfold < ../data/" + inputfile_files + ".random.fasta > ../data/" + inputfile_files + ".random.fold"
    os.system(cmd)

    cmd = "perl miPred/genRNARandomStats.pl -n 1000 -i ../data/" + inputfile_files + ".random.fold -o ../data/" + inputfile_files + ".zdata -m ../data/" + inputfile_files + ".fold"
    os.system(cmd)

    #----------------------------------------------------------------------------------

    #calculating MFEI1, MFEI2, MFEI3, MFEI4
    cmd = "java mfe14 ../data/" + inputfile_files
    os.system(cmd)

    cmd = "java mfe23 ../data/" + inputfile_files
    os.system(cmd)

    #----------------------------------------------------------------------------------
    
    #calculating basepair-related features
    cmd = "java bpcount1 ../data/" + inputfile_files
    os.system(cmd)
    
    cmd = "java bpcount2 ../data/" + inputfile_files
    os.system(cmd)

    #----------------------------------------------------------------------------------

    #calculating RNAfold-related features
    cmd = "./RNAfold/RNAfold -p2 < ../data/" + inputfile_files + " > ../data/" + inputfile_files + ".RNAfold1"
    os.system(cmd)
    
    cmd = "java RNAfoldfilter ../data/" + inputfile_files
    os.system(cmd)
    
    cmd = "./hybrid-ss-min ../data/" + inputfile_files + " -o ../data/" + inputfile_files
    os.system(cmd)

    #----------------------------------------------------------------------------------

    #calculating mfold-related features
    cmd = "perl melt.pl ../data/" + inputfile_files + " > ../data/" + inputfile_files + ".mfold"
    os.system(cmd)

    cmd = "java Mfoldfilter ../data/" + inputfile_files
    os.system(cmd)
 
    #----------------------------------------------------------------------------------
    
    #filtering all features
    cmd = "java filter_all " + inputfile_files
    os.system(cmd)
    
    cmd = "java format1 " + inputfile_files
    os.system(cmd)
    
    #----------------------------------------------------------------------------------

    #copying features file to ../results/ + inputfilename
    cmd = "mv ../data/selected." + inputfile_files + ".-21.features ../results/" + inputfile + "/selected." + str(count) + ".-21.features"
    os.system(cmd) 

    try:
      f3.close()
    except:
      pass
    
    input = "../results/" + inputfile + "/selected." + str(count) + ".-21.features"
    
    #----------------------------------------------------------------------------------

    #checking if features are calculated completely and writing to 'selected_21_micropred_features'
    f_input = open(input, 'r')
    try:
      file = f_input.readlines()[0]
    except:
      file = '?  ' * 20
      file += '?\n'
    if len(file.split('  ')) == 21:
      file = file.replace('  ', ',')
      out_main.write(file) 
    else:
      write = '?,' * 20
      write += '?\n'
      out_main.write(write)
    f_input.close()

    cmd = 'rm ' + input
    os.system(cmd)

#----------------------------------------------------------------------------------

cmd = "rm ../data/" + inputfile_files + "*"
os.system(cmd)
cmd = 'rm ../data/*.features'
os.system(cmd)
cmd = 'rm *.ps'
os.system(cmd)


  
