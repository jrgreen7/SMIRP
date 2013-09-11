import sys, os

if len(sys.argv) != 3:
	print "Incorrect number of parameters provided!"
	print "USAGE:\nclassify modelFile featuresFile\n\nPARAMETERS:"
	print "modelFile - file containing classification model. It can be a model trained by the user in Weka software or one of the precomputed models:"
	print "\thuman.model\n\tarabidopsis.model\n\tanimal.model\n\tplant.model\n\tvirus.model"
	print	"featuresFile - arff file with representation of sequences to be classified (21+7 features).\n"

else:	
	modelFile = "classifier/" + sys.argv[1]
	inputFile = "results/" + sys.argv[2]
	outFile =  inputFile + ".predictions"
	
	cmd = "java -classpath classifier" + os.pathsep + "classifier/weka.jar" + os.pathsep + "classifier/ROCSelect.jar huntmi.HuntMi " + modelFile + " " + inputFile + " " + outFile
	os.system(cmd)