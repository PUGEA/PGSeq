#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  


import os
import time
import calculateModelData
import gc



def generateData(AlignmentNames,InputName,targetGeneFile,ReadType,ReadLen):

    currentPath = os.getcwd()

    annotationPath = os.path.join(currentPath,'Annotation')
    AnnoGeneLabel = os.path.join(annotationPath,InputName+'.Gene.Label.Isoform')
    AnnoExonSplit = os.path.join(annotationPath,InputName+'.Exon.Split')
    AnnoGeneMap = os.path.join(annotationPath,InputName+'.Gene.Map.Isoform')
    AnnoExonLen = os.path.join(annotationPath,InputName+'.Exon.Len')


    tempDataPath = os.path.join(currentPath,'Model.Temp')
    modelDataPath = os.path.join(currentPath,'Model.Data')
    inputDataFiles = []
    AlignmentReadNo = []

    ReadLen = int(ReadLen)


    for index in xrange(len(AlignmentNames)):
        alignName = AlignmentNames[index]

        print '        ',index+1,'.', alignName

        inputAlignPath = os.path.join(tempDataPath,alignName)

        # print 'step CalculateProbability....'

        inputPath = os.path.join(inputAlignPath,'Data.Covert')
        inputFile = os.path.join(inputPath,alignName+'.covert')
        outputPath = os.path.join(inputAlignPath,'Data.Prob')
        outputFile = os.path.join(outputPath,alignName)
        ReadNo = calculateModelData.CalculateProbability(inputFile,outputFile,ReadType)
        AlignmentReadNo.append(ReadNo)



        # print 'step ExtractGeneData...'


        inputFile = outputFile
        outputPath = os.path.join(inputAlignPath,'Gene.LocRel')
        calculateModelData.ExtractGeneData(inputFile,outputPath,targetGeneFile,AnnoGeneLabel,ReadType)

        # print 'step CalculateLocation...'


        inputPath = outputPath
        outputPath = os.path.join(inputAlignPath,'Gene.LocAbs')
        inputDataFiles.append(outputPath)
        calculateModelData.CalculateLocation(inputPath,outputPath,targetGeneFile,ReadLen,ReadType,AnnoExonSplit,AnnoExonLen)   

        gc.collect()

    print '   3.2. Obtaining the Model Data .... '
 
    # print 'statisticGeneMap...' 

    outputPath = os.path.join(modelDataPath,'Gene.Map')
    calculateModelData.statisticGeneMap(targetGeneFile,outputPath,AnnoGeneMap)

    # print 'statisticGeneLength...' 

    outputPath = os.path.join(modelDataPath,'Gene.Length')
    calculateModelData.statisticGeneLength(targetGeneFile,outputPath,AnnoExonLen)

    # print 'statisticGeneData...' 

    outputPath = os.path.join(modelDataPath,'Gene.Data')
    calculateModelData.statisticGeneData(targetGeneFile,inputDataFiles,outputPath,AnnoExonLen)

    f = open(os.path.join(modelDataPath,'Alignment.Read.No'),'w')
    for i in range(len(AlignmentReadNo)):
        f.write(str(AlignmentReadNo[i])+'\t')
    f.close()







