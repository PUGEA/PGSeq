#! /usr/bin/env python 
from __future__ import division  


import os 
import sys, getopt
import example
import multiprocessing as mp
import time
import writeExpression
import createAlignmentFolder
import analyseAnnotation

def usage():
  print "Usage:"
  print "    calculateExpression -a AnnotationPrefix  -s SelectedGenes "
  print "options:"
  print '   -a/--AnnotationPrefix <text>   the abbreviative name is same as the outputName of preprocessAnnotation, eg:ensGene'
  print '   -s/--SelectedGenes    <path>   input objective gene'
  print '   -l/--log                       the natural logarithmic expression to detecting different expressed genes or isoforms'              
  print '   -h/--help                      display this help and exit'
  print '   -v/--version                   output version information and exit\n'
  sys.exit(1)

def error():
    usage()



try:
    opts, args = getopt.getopt(sys.argv[1:], "hvla:s:", ["help","version", "AnnotationPrefix","SelectedGenes","log"])
except getopt.GetoptError:
    error()
    sys.exit()


targetGene = 'default'
input_name = 'default'
log_flag =0

for op, value in opts:
    if op in ["-a","--AnnotationPrefix"]:
        input_name = value
    elif op in ['-s','--SelectedGenes']:
        targetGene = value
    elif op in ['-l','--log']:
        log_flag = 1
    elif op in ["-v","--version"]:
        print 'The version of PBSeq is v0.2_alpha.'
        sys.exit()
    elif op in ["-h","--help"]:
        usage()        
        sys.exit()

if targetGene == 'default':                            # now, the targetGene will set 'default'
    currentPath = os.getcwd()
    targetDir = os.path.join(currentPath,'Annotation')
    targetGene = os.path.join(targetDir,input_name)
    targetGeneFile = targetGene+'.Gene.Info'  
else:
    targetDir = os.getcwd()
    outputDir = os.path.join(targetDir,'Annotation')
    targetGeneFile = analyseAnnotation.SelectedGenesInfo(targetGene,input_name,outputDir)

parametersList = [input_name,targetGeneFile]

if 'default' in parametersList:
    usage()
    sys.exit()
else:
    pass

f= open(targetGeneFile,'r')
geneNameList = []
dictIsoNo = {}
dictTransName = {}
for line in f:
    line = line.rstrip()
    line = line.split('\t')
    geneName = line[0]
    isoNo = int(line[1])
    transName = line[3:]
    geneNameList.append(geneName)
    dictIsoNo[geneName]=isoNo
    dictTransName[geneName]=transName
f.close()




currentPath = os.getcwd()
modelDataPath = os.path.join(currentPath,'Model.Data')
resultDataPath = os.path.join(currentPath,'PGSeq.Results')

f = open(os.path.join(modelDataPath,'Alignment.Read.No'),'r')
Data = f.readline()
f.close()
Data = Data.rstrip()
AlignmentReadNo = Data.split('\t')
repNo = len(AlignmentReadNo)
NormalizeConst = [0 for i in range(repNo)]
for i in range(repNo):
    NormalizeConst[i] = int(AlignmentReadNo[i])/int(AlignmentReadNo[0])
example.getNormc(NormalizeConst)


AnnoDir = os.path.join(currentPath,'Annotation')
AnnoExonFile = os.path.join(AnnoDir,input_name+'.Exon.Len')
dictExonNo = {}
f = open(AnnoExonFile,'r')
for line in f:
    line = line.rstrip()
    line = line.split('\t')
    geneName = line[0]
    exonNo = line[1]
    dictExonNo[geneName] = int(exonNo)
f.close()




def geneCalculateExpression(geneIndex):
    geneName = geneNameList[geneIndex]
    isoNo = dictIsoNo[geneName]
    exonNo = dictExonNo[geneName]

    geneDataPath = os.path.join(modelDataPath,'Gene.Data')
    geneDataFile = os.path.join(geneDataPath,geneName)

    geneMapPath = os.path.join(modelDataPath,'Gene.Map')
    geneMapFile = os.path.join(geneMapPath,geneName)

    geneLenPath = os.path.join(modelDataPath,'Gene.Length')
    geneLenFile = os.path.join(geneLenPath,geneName)

    example.getgenedata(geneDataFile,geneMapFile,geneLenFile,isoNo,exonNo)   
    expression =example.calexpression(isoNo,repNo)
    dictExpression[geneName]=expression  


print ''
# print 'calculateExpression....',time.ctime()
print '1. Calculating Expression .... '
if __name__ == "__main__":

    manager = mp.Manager()
    dictExpression = manager.dict()
    pool = mp.Pool(processes=mp.cpu_count()*2)  
    for index in xrange(len(geneNameList)):
      pool.apply_async(geneCalculateExpression, (index, ))
    pool.close()
    pool.join()

    print '2. Writing Expression ....'

    writeExpression.writeExpression(targetGeneFile,dictExpression,repNo,log_flag,input_name)


createAlignmentFolder.deleteDataFolder()

print 'Done.'



















