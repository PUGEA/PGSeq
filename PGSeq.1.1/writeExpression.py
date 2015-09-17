#! /usr/bin/env python

import os
import gc
import time 
import numpy as np
import multiprocessing as mp

SamplingNo = 20000
const = np.log(2) 

def covertSequence(Seq):
    Seq = Seq.strip()
    Seq = Seq.split(',')
    Seq.pop()
    l = len(Seq)
    SeqCovert = []
    for i in range(l):
        SeqCovert.append(int(Seq[i]))
    return np.array(SeqCovert)


def normGeneSample(geneName,geneMean,geneStd,repNo,index):
  geneLogMean = []
  geneLogStd = []

  for i in xrange(repNo):
      np.random.seed(index)
      samples = np.random.normal(geneMean[i],geneStd[i],SamplingNo)
      samplesPos = samples[samples>0]
      samplesPosLog = np.log(samplesPos)
      geneLogMean.append(np.mean(samplesPosLog))
      geneLogStd.append(np.std(samplesPosLog))
  dictGeneLogMean[geneName]=np.array(geneLogMean)
  dictGeneLogStd[geneName]=np.array(geneLogStd)

def normTransSample(transName,transMean,transStd,repNo,index):
  transLogMean = []
  transLogStd = []

  for i in xrange(repNo):
      np.random.seed(index)
      samples = np.random.normal(transMean[i],transStd[i],SamplingNo)
      samplesPos = samples[samples>0]
      samplesPosLog = np.log(samplesPos)
      transLogMean.append(np.mean(samplesPosLog))
      transLogStd.append(np.std(samplesPosLog))    
  dictTransLogMean[transName]=np.array(transLogMean)
  dictTransLogStd[transName]=np.array(transLogStd)


def writeExpression(targetGeneFile,dictExpression,repNo,logFlag,input_name):

    currentPath = os.getcwd()
    resultDataPath = os.path.join(currentPath,'PGSeq.Results') 

    f= open(targetGeneFile,'r')
    geneNameList = []
    dictIsoNo = {}
    dictTransName = {}
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        geneName = line[0]
        isoNo = int(line[1])
        transNameList = line[3:]

        geneNameList.append(geneName)
        dictIsoNo[geneName]=isoNo
        dictTransName[geneName]=transNameList
    f.close()
    geneNo = len(geneNameList)

    dictTransMean={}
    dictTransVar={}
    dictGeneMean={}
    dictGeneVar={}

    for index in xrange(geneNo):
        geneName = geneNameList[index]
        isoNo = dictIsoNo[geneName]
        transNameList = dictTransName[geneName]

        alphaNo = isoNo*repNo
        geneExpress = dictExpression[geneName]

        geneMean = np.array(geneExpress[2*alphaNo:2*alphaNo+repNo])
        geneVar = np.array(geneExpress[2*alphaNo+repNo:])
        dictGeneMean[geneName]= geneMean
        dictGeneVar[geneName]= geneVar

        transMean = np.array(geneExpress[0:alphaNo])
        transVar = np.array(geneExpress[alphaNo:2*alphaNo])
        transMean.shape = (repNo,isoNo)
        transVar.shape = (repNo,isoNo)

        for i in xrange(isoNo):
            transName = dictTransName[geneName][i]
            dictTransMean[transName] = transMean[:,i]
            dictTransVar[transName] = transVar[:,i]

    del dictExpression
    gc.collect()

    f_transMean=open(os.path.join(resultDataPath,'isoforms.mean'),'w')
    f_transVar=open(os.path.join(resultDataPath,'isoforms.std'),'w')
    f_geneMean=open(os.path.join(resultDataPath,'genes.mean'),'w')
    f_geneVar=open(os.path.join(resultDataPath,'genes.std'),'w')

    for index in xrange(geneNo):
        geneName = geneNameList[index]
        isoNo = dictIsoNo[geneName]
        transNameList = dictTransName[geneName]

        geneMean = dictGeneMean[geneName]
        geneVar = dictGeneVar[geneName]  

        f_geneMean.write(geneName+'\t')
        f_geneVar.write(geneName+'\t')
        for i in xrange(repNo):    
            f_geneMean.write(str(round(geneMean[i]/const,6))+'\t')
            f_geneVar.write(str(round(np.sqrt(geneVar[i])/const,6))+'\t')
        f_geneMean.write('\n')
        f_geneVar.write('\n')


        for i in xrange(isoNo): 
            transName = transNameList[i]
            transMean = dictTransMean[transName]
            transVar = dictTransVar[transName]

            f_transMean.write(geneName+'\t')
            f_transVar.write(geneName+'\t')
            f_transMean.write(transName+'\t')
            f_transVar.write(transName+'\t')        
            for j in xrange(repNo):
                f_transMean.write(str(round(transMean[j]/const,6))+'\t')
                f_transVar.write(str(round(np.sqrt(transVar[j])/const,6))+'\t')
            f_transMean.write('\n')
            f_transVar.write('\n')

    f_transMean.close()
    f_transVar.close()
    f_geneMean.close()
    f_geneVar.close()



    if logFlag == 1:
        # print 'Sampling...', time.ctime()

        # AnnoDir = os.path.join(currentPath,'Annotation')
        # AnnoExonFile = os.path.join(AnnoDir,input_name+'.Exon.Len')  
        # dictExonLen = {}
        # dictGeneLen = {}        
        # f = open(AnnoExonFile,'r')
        # for line in f:
        #     line = line.rstrip()
        #     line = line.split('\t')
        #     geneName = line[0]
        #     exonLen = covertSequence(line[2])
        #     dictExonLen[geneName] = exonLen
        #     dictGeneLen[geneName] = sum(exonLen)
        # f.close()   

        # AnnoTransMap = os.path.join(AnnoDir,input_name+'.Gene.Map.Isoform')
        # dictTransLen = {}  
        # f = open(AnnoTransMap,'r')
        # for line in f:
        #     line = line.rstrip()
        #     line = line.split('\t')
        #     geneName = line[0]
        #     transName =line[1]
        #     transMap = covertSequence(line[2])
        #     exonLen = dictExonLen[geneName]
        #     dictTransLen[transName] = sum(transMap*exonLen)
        # f.close()

        global dictGeneLogMean
        global dictGeneLogStd
        global dictTransLogMean
        global dictTransLogStd

        manager = mp.Manager()        
        dictGeneLogMean= manager.dict()
        dictGeneLogStd = manager.dict()
        dictTransLogMean = manager.dict()
        dictTransLogStd = manager.dict()  

        pool = mp.Pool(processes=mp.cpu_count()*2)  
        # pool = mp.Pool(processes=3)
        for index in xrange(geneNo):
            geneName = geneNameList[index]
            geneMean = dictGeneMean[geneName]
            geneStd = np.sqrt(dictGeneVar[geneName])
            pool.apply_async(normGeneSample, (geneName,geneMean,geneStd,repNo,index,))
        pool.close()
        pool.join()   

        pool = mp.Pool(processes=mp.cpu_count()*2)  
        # pool = mp.Pool(processes=3)
        for index in xrange(geneNo):
            geneName = geneNameList[index]
            isoNo = dictIsoNo[geneName]
            transNameList = dictTransName[geneName]
            for i in xrange(isoNo):
              transName = transNameList[i]              
              transMean = dictTransMean[transName]
              transStd = np.sqrt(dictTransVar[transName])
              pool.apply_async(normTransSample, (transName,transMean,transStd,repNo,index,))
        pool.close()
        pool.join()   


        f_translogMean=open(os.path.join(resultDataPath,'isoforms.mean.log'),'w')
        f_translogVar=open(os.path.join(resultDataPath,'isoforms.std.log'),'w')
        f_genelogMean=open(os.path.join(resultDataPath,'genes.mean.log'),'w')
        f_genelogVar=open(os.path.join(resultDataPath,'genes.std.log'),'w')

        for index in xrange(geneNo):
            geneName = geneNameList[index]
            isoNo = dictIsoNo[geneName]
            transNameList = dictTransName[geneName]
            geneLogMean = dictGeneLogMean[geneName]
            geneLogStd = dictGeneLogStd[geneName] 


            f_genelogMean.write(geneName+'\t')
            f_genelogVar.write(geneName+'\t') 
            for i in xrange(int(repNo)):           
                f_genelogMean.write(str(round(geneLogMean[i]/const,6))+'\t')
                f_genelogVar.write(str(round(geneLogStd[i]/const,6))+'\t')        
            f_genelogMean.write('\n')
            f_genelogVar.write('\n')


            for i in xrange(isoNo):
              transName = transNameList[i]
              transLogMean = dictTransLogMean[transName]
              transLogStd = dictTransLogStd[transName]

              f_translogMean.write(geneName+'\t')
              f_translogVar.write(geneName+'\t')
              f_translogMean.write(transName+'\t')
              f_translogVar.write(transName+'\t')
              for j in xrange(int(repNo)):  
                  f_translogMean.write(str(round(transLogMean[j]/const,6))+'\t')
                  f_translogVar.write(str(round(transLogStd[j]/const,6))+'\t')
              f_translogMean.write('\n')
              f_translogVar.write('\n') 

        f_translogMean.close()
        f_translogVar.close()
        f_genelogMean.close()
        f_genelogVar.close() 

    else:
      pass


    # print 'Done.', time.ctime()

