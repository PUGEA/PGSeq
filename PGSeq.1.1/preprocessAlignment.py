 #!/usr/bin/env python


import os, os.path
import sys, getopt
import analyseAnnotation
import createAlignmentFolder
import covertAlignment
import generateData


def usage():
  print "Usage:"
  print "    preprocessAlignment  -t AnnotationType -a AnnotationFile -o AnnotationPrefix -d AlignmentFiles -s SelectedGenes"
  print "options:"  	
  print "   -t/--AnnotationType   <int>    supported four annotation types: refGene: 1, ensGene: 2, knownGene: 3 and Ensembl: 4"
  print '   -a/--AnnotationFile   <path>   input the reference annotation file, eg: ensGene.txt'
  print '   -o/--AnnotationPrefix <text>   the abbreviative name of output file, eg: ensGene'
  print '   -d/--AlignmentFiles   <path>   input one or a list of alignment files from Bowtie2, eg: data1.sam'
  print '   -s/--SelectedGenes    <paht>   input objective gene'
  print '   -h/--help                      display this help and exit'
  print '   -v/--version                   output version information and exit\n'

  print "Example:"
  print "    Supported four types: refGene, ensGene, knownGene and Ensembl"
  print "    preprocessAlignment.py -t 1 -a refGene.txt  -o refGene -d data1.sam,data2,sam"
  print "    preprocessAlignment.py -t 2 -a ensGene.txt  -o ensGene -d data1.sam,data2,sam"
  print "    preprocessAlignment.py -t 3 -a knownGene.txt,knownIsoforms.txt  -o knownGene -d data1.sam,data2,sam"
  print "    preprocessAlignment.py -t 4 -a Homo_sapiens.GRCh37.67.gtf  -o Ensembl -d data1.sam,data2,sam"
  sys.exit(1)

  	

def error():
    usage()

try:
    opts, args = getopt.getopt(sys.argv[1:], "hvt:a:o:d:s:", ["help","version","AnnotationType","AnnotationFile","AnnotationPrefix","AlignmentFiles","SelectedGenes"])
except getopt.GetoptError:
    usage()
    sys.exit()



          
targetGene = 'default'
prefix_name = 'default'
alignment_paths = 'default'

for op, value in opts:
    if op in ["-a","--AnnotationFile"]:
        annotation_file = value
    elif op in ["-t","--AnnotationType"]:
        typeflag = value
    elif op in ["-o","--AnnotationPrefix"]:
        prefix_name = value
    elif op in ["-d","--AlignmentFiles"]:
        alignment_paths = value        
    elif op in ['-s','--SelectedGenes']:
        targetGene = value        
    elif op in ["-v","--version"]:
        print 'The version of PBSeq is v0.2_alpha.'
        sys.exit()
    elif op in ["-h","--help"]:
        usage()        
        sys.exit()

if len(opts)<4:
    usage()
    sys.exit()
else:
    pass

AnnoType = {'1':'refGene','2':'ensGene','3':'knownGene','4':'Ensembl'}
annotation_type = AnnoType[typeflag]



parametersList = [alignment_paths,prefix_name]

if 'default' in parametersList:
    error()
    usage()
    sys.exit()
else:
    pass


alignment_paths = alignment_paths.rstrip()
alignment_paths = alignment_paths.split(',')
Alignment_Names = []
for i in xrange(len(alignment_paths)):
  	alignment_path = alignment_paths[i]
  	alignment_path = alignment_path.split('/')
  	alignfile = alignment_path[-1]

  	alignfile = alignfile.split('.')
  	alignfile = alignfile[0:len(alignfile)-1]
  	alignName = alignfile[0]
  	for j in xrange(1,len(alignfile)):
    		alignName = alignName+'.'+alignfile[j]
  	Alignment_Names.append(alignName)
# print Alignment_Names

# print alignment_paths[0]
f = open(alignment_paths[0],'r')
for line in f:
	line = line.rstrip()
	line = line.split('\t')
	if line[2]=='*':
		pass
	else:
		if line[6] =='=':
			ReadType = 'Paired'
		else:
			ReadType = 'Single'
		ReadLen = len(line[9])
		break
f.close()

# print ReadType, ReadLen



print ''
print '1. Creating Folders .... '
createAlignmentFolder.createFolder(Alignment_Names)


targetDir = os.getcwd()
targetFolder = os.path.join(targetDir,'Annotation')
if os.path.exists(targetFolder):
    pass
else:
    os.mkdir(targetFolder)
outputDir = targetFolder

print "2. Preprocessing Annotation Files .... "
analyseAnnotation.CovertAnnotation(annotation_type,annotation_file,prefix_name,outputDir)
analyseAnnotation.AnalyseAnnotation(annotation_file,prefix_name,outputDir)


if targetGene == 'default':                            # now, the targetGene will set 'default'
    currentDir = os.getcwd()
    targetDir = os.path.join(currentDir,'Annotation')
    targetGene = os.path.join(targetDir,prefix_name)
    targetGeneFile = targetGene+'.Gene.Info'  
else:
    targetGeneFile = analyseAnnotation.SelectedGenesInfo(targetGene,prefix_name,outputDir)
     



print '3. Preprocessing Alignment Files .... '
print '   3.1. Coverting Alignment Files .... '
covertAlignment.covertAlignment(annotation_type,alignment_paths,Alignment_Names,prefix_name,ReadType)


print '   3.2. Preprocessing Alignment Files .... '

generateData.generateData(Alignment_Names,prefix_name,targetGeneFile,ReadType,ReadLen)



createAlignmentFolder.deleteTempFolder()

print 'Done.'