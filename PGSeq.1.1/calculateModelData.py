#! /usr/bin/env python

import pp
import os
import random as rd
import numpy as np
import gc


def CalculateProbability(InputFile,OutputFile,ReadType):

    DictRead = {}

    f = open(InputFile,'r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        read = line[0]
        if ReadType == 'Single':
            gene = line[3]
        if ReadType == 'Paired':
            gene = line[4] 

        if read not in DictRead:
            genelist = [gene]
            DictRead[read] = genelist
        else:
            # genelist = DictRead[read]
            if gene in genelist:
                pass
            else:
                genelist.append(gene)
                DictRead[read] = genelist
    f.close()

    GeneUniqRead = {}

    for key in DictRead.iterkeys():
        read = key
        genelist = DictRead[key]

        if len(genelist) == 1:
            gene = genelist[0]
            if gene not in GeneUniqRead:
               GeneUniqRead[gene] = 1
            else:
               GeneUniqRead[gene] = GeneUniqRead[gene] + 1
        else:
            pass


    for key in DictRead.iterkeys():
        read = key
        genelist = DictRead[key]
        geneno = len(genelist)

        if geneno == 1:
            geneRate = [1.0]
            DictRead[read] = [genelist,geneRate]
        else:
            countlist = []
            for i in xrange(geneno):
                gene = genelist[i]
                if gene in GeneUniqRead:
                    count = GeneUniqRead[gene]
                else:
                    count = 0 
                countlist.append(count)

            if sum(countlist) ==0:                   # 
                rd.seed(123456)
                rdno = rd.randint(0,geneno-1)
                geneRate = [0.0 for j in xrange(geneno)]
                geneRate[rdno] = 1.0
            else:
                MaxGene = max(countlist)
                MaxInd = countlist.index(MaxGene)
                geneRate = [0.0 for j in xrange(geneno)]
                geneRate[MaxInd] = 1.0

            DictRead[read] = [genelist,geneRate]          

    

    f = open(InputFile,'r')
    fo = open(OutputFile+'.Prob','w')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        
        if ReadType == 'Single':
            read, loc, trans, gene = line
        if ReadType == 'Paired':
            read, loc1, loc2, trans, gene = line



        genelist, geneRate = DictRead[read]
        geneInd = genelist.index(gene)
        readRate = geneRate[geneInd]
        
        if ReadType == 'Single':
            fo.write(read+'\t'+gene+'\t'+trans+'\t'+loc+'\t'+str(readRate)+'\n')
        if ReadType == 'Paired':
            fo.write(read+'\t'+gene+'\t'+trans+'\t'+loc1+'\t'+str(loc2)+'\t'+str(readRate)+'\n')       
    
    f.close()
    fo.close()

    flab = open(OutputFile+'.Read','w')
    readNo = 0
    for read in DictRead.iterkeys():
        flab.write(read+'\t'+str(readNo)+'\n')
        readNo = readNo+1
    flab.close()

    return readNo


def ExtractGeneData(InputFile,OutputPath,ObjectGene,FileLabel,ReadType):


    DictTransLab = {}
    DictGeneLab = {}
    f = open(FileLabel,'r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        gene,lab, trans = line

        if gene not in DictGeneLab:
            DictGeneLab[gene] = [trans]
        else:
            temp  = DictGeneLab[gene]
            temp.append(trans)
            DictGeneLab[gene] =temp        
        DictTransLab[trans] = int(lab)
    f.close()

    DictReadLab = {}
    DictReadLabC = {}
    f = open(InputFile+'.Read','r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        read, lab = line     
        DictReadLab[read] = lab
        DictReadLabC[lab] = read
    f.close()

    DictGeneData = {}    
    if ReadType == 'Single':   
        f = open(InputFile+'.Prob','r') 
        for line in f:
            line = line.rstrip()
            line = line.split('\t')
            read,gene,trans,loc,prob = line

            readlab = DictReadLab[read]
            translab = DictTransLab[trans]

           
            if gene not in DictGeneData:
                readlist = [readlab]
                problist = [prob]
                translist = [translab]
                loclist = [loc]
                DictGeneData[gene]= readlist,problist,translist,loclist
            else:
                readlist,problist,translist,loclist = DictGeneData[gene]
                readlist.append(readlab)
                problist.append(prob)
                loclist.append(loc)
                translist.append(translab)                
                DictGeneData[gene]= readlist,problist,translist,loclist
        f.close()
    elif ReadType == 'Paired':
        f = open(InputFile+'.Prob','r')         
        for line in f:
            line = line.rstrip()
            line = line.split('\t')
            read,gene,trans,loc1,loc2,prob = line
    
            readlab = DictReadLab[read]
            translab = DictTransLab[trans]
            
            if gene not in DictGeneData:
                readlist = [readlab]
                problist = [prob]
                translist = [translab]
                loc1list = [loc1]
                loc2list = [loc2]
                DictGeneData[gene]= readlist,problist,translist,loc1list,loc2list
            else:
                readlist,problist,translist,loc1list,loc2list = DictGeneData[gene]
                readlist.append(readlab)
                problist.append(prob)
                loc1list.append(loc1)
                loc2list.append(loc2)
                translist.append(translab)                
                DictGeneData[gene]= readlist,problist,translist,loc1list,loc2list
        f.close()
    else:
        pass
                


    if ReadType == 'Single':
        f = open(ObjectGene,'r')
        for line in f:
            line = line.rstrip()
            line = line.split('\t')
            gene = line[0].rstrip()           

            if gene in DictGeneData:
                readlist,problist,translist,loclist = DictGeneData[gene]
                transName = DictGeneLab[gene]
    
                fo = open(OutputPath+'/'+gene,'w')
                for i in range(len(readlist)):
                    read = DictReadLabC[readlist[i]]
                    trans = transName[translist[i]-1]
                    fo.write(read+'\t'+trans+'\t'+loclist[i]+'\t'+problist[i]+'\n')
                fo.close()
            else:
                pass
        f.close()
            
    elif ReadType == 'Paired':
        f = open(ObjectGene,'r')
        for line in f:
            line = line.rstrip()
            line = line.split('\t')
            gene = line[0].rstrip()  

            if gene in DictGeneData:
                readlist,problist,translist,loc1list,loc2list = DictGeneData[gene]
                transName = DictGeneLab[gene]
    
                fo = open(OutputPath+'/'+gene,'w')
                for i in range(len(readlist)):
                    read = DictReadLabC[readlist[i]]
                    trans = transName[translist[i]-1]
                    fo.write(read+'\t'+trans+'\t'+loc1list[i]+'\t'+loc2list[i]+'\t'+problist[i]+'\n')
                fo.close()
            else:
                pass   
        f.close()               
    else:
        pass

    # os.remove(InputFile+'.Read')
    # os.remove(InputFile+'.Prob')

    del DictReadLab
    gc.collect()
    del DictReadLabC
    gc.collect()   
    del DictGeneData
    gc.collect()





def Covert_Sequence(Seq):
    l = len(Seq)
    SeqCovert = []
    for i in range(l):
        SeqCovert.append(int(Seq[i]))
    return SeqCovert



    
def Calculate_GeneLocation(GeneName,IsoNo,InputFile,OutputFile,ReadLen,Dict_Gene_Info,Dict_Iso_Info,Dict_Exon_Info,ReadType):    
    
    f_in = open(InputFile,'r')
    f_out = open(OutputFile,'w')
    
    GeneExonNo,GeneExonLab,GeneExonLen = Dict_Gene_Info[GeneName]
    for line in f_in:
        line = line.rstrip()
        line = line.split('\t')
        ReadId = line[0]
        IsoName = line[1]
        if ReadType == 'Single':
            Loc1 = int(line[2])
            Loc2 = Loc1 + ReadLen
            Weight = line[3]
        elif ReadType == 'Paired':                
            Loc1 = int(line[2])
            Loc2 = int(line[3])+ReadLen
            Weight = line[4]  
        else:
            pass

            
        
        IsoExonNo,IsoExonLab,IsoExonLen = Dict_Iso_Info[IsoName]
    
        
        if Loc1 > IsoExonLen[-1] or Loc2+ReadLen > IsoExonLen[-1]:
            pass    
        else:   
            if IsoExonNo == 1:
                Rel_Ind1 = IsoExonLab[0]
                Rel_Ind2 = IsoExonLab[0]            
                Rel_Loc1 = Loc1
                Rel_Loc2 = Loc2
                if Rel_Ind1 == 1:
                    Abs_Loc1 = Rel_Loc1
                    Abs_Loc2 = Rel_Loc2
                else:
                    Abs_Loc1 = GeneExonLen[Rel_Ind1-2]+Rel_Loc1
                    Abs_Loc2 = GeneExonLen[Rel_Ind2-2]+Rel_Loc2 
                
            else:
                for i in range(IsoExonNo-1):
                    if Loc1 <= IsoExonLen[0]:
                        Rel_Ind1 = IsoExonLab[0]
                        Rel_Loc1 = Loc1
                        if Rel_Ind1 ==1:
                            Abs_Loc1 = Rel_Loc1
                        else:
                            Abs_Loc1 = GeneExonLen[Rel_Ind1-2]+Rel_Loc1
                            
                        if Loc2<= IsoExonLen[0]:
                            Rel_Ind2 = IsoExonLab[0]
                            Rel_Loc2 = Loc2
                            if Rel_Ind2 ==1:
                                Abs_Loc2 = Rel_Loc2
                            else:
                                Abs_Loc2 = GeneExonLen[Rel_Ind2-2]+Rel_Loc2  
                        
                        else:
                            for i in range(IsoExonNo-1):
                                if Loc2 <= IsoExonLen[i+1] and Loc2 > IsoExonLen[i]:
                                    Rel_Ind2 = IsoExonLab[i+1]
                                    Rel_Loc2 = Loc2 - IsoExonLen[i]
                                    Abs_Loc2 = GeneExonLen[Rel_Ind2-2]+Rel_Loc2
                                    break
                                else:
                                    continue
                    elif Loc1 <= IsoExonLen[i+1] and Loc1 > IsoExonLen[i]:
                        Rel_Ind1 = IsoExonLab[i+1]
                        Rel_Loc1 = Loc1 - IsoExonLen[i]
                        Abs_Loc1 = GeneExonLen[Rel_Ind1-2]+Rel_Loc1
                        
                        for j in range(i,IsoExonNo-1):
                            if Loc2 <= IsoExonLen[j+1] and Loc2 > IsoExonLen[j]:
                                Rel_Ind2 = IsoExonLab[j+1]
                                Rel_Loc2 = Loc2 - IsoExonLen[j]
                                Abs_Loc2 = GeneExonLen[Rel_Ind2-2]+Rel_Loc2
                                break
                            else:
                                continue                    
                        break                    
                    else:
                        continue
                    
                
            
            f_out.write(ReadId+'\t'+IsoName+'\t'+str(Weight)+'\t'+str(Abs_Loc1)+'\t'+str(Abs_Loc2)+'\t')
            FragCrossExon1 = IsoExonLab.index(Rel_Ind1)
            FragCrossExon2 = IsoExonLab.index(Rel_Ind2)
            FragCrossExon = IsoExonLab[FragCrossExon1:FragCrossExon2+1]
            for k in range(len(FragCrossExon)):
                f_out.write(str(FragCrossExon[k])+',')
            f_out.write('\t')
            
            ExonLenList = Dict_Exon_Info[GeneName]
            CrossExonNo = len(FragCrossExon)
            FragExonLen = []

            if CrossExonNo==1:
                FragExonLen.append(Loc2-Loc1)
            elif CrossExonNo == 2:
                FragExonLen.append(GeneExonLen[FragCrossExon[0]-1]-Abs_Loc1)
                FragExonLen.append(Abs_Loc2-GeneExonLen[FragCrossExon[1]-2])
            else:
                FragExonLen.append(GeneExonLen[FragCrossExon[0]-1]-Abs_Loc1)
                
                for k in FragCrossExon[1:CrossExonNo-1]:
                    FragExonLen.append(ExonLenList[k-1])
                FragExonLen.append(Abs_Loc2-GeneExonLen[FragCrossExon[-1]-2])
                
            for k in range(len(FragCrossExon)):
                f_out.write(str(FragExonLen[k])+',')
            f_out.write('\n')        
    f_in.close()
    f_out.close()



def CalculateLocation(InputPath,OutputPath,ObjectGene,ReadLen,ReadType,FileExon,FileLen):

    f_ref = open(FileExon,'r')
    Dat = [line.strip() for line in f_ref.readlines()];
    Dat = [line.split('\t') for line in Dat];
    f_ref.close()

    Dict_Gene_Info= {}
    Dict_Iso_Info = {}
    
    i=0
    while i < len(Dat):
        Gene = Dat[i][0]
        IsoNo = int(Dat[i][1])
        GeneExonNo = int(Dat[i][2])
        GeneExonLab = Dat[i][3].split(',')
        GeneExonLab.pop()
        GeneExonLab = Covert_Sequence(GeneExonLab)
        GeneExonLen = Dat[i][4].split(',')
        GeneExonLen.pop()
        GeneExonLen = Covert_Sequence(GeneExonLen)
        
        Dict_Gene_Info[Gene] = GeneExonNo,GeneExonLab,GeneExonLen
        
        IsoDat = Dat[i+1:i+IsoNo+1]
        for j in xrange(IsoNo):
            GeneName = IsoDat[j][0]
            IsoName = IsoDat[j][1]
            IsoExonNo = int(IsoDat[j][2])
            IsoExonLab = IsoDat[j][3].split(',')
            IsoExonLab.pop()
            IsoExonLab = Covert_Sequence(IsoExonLab)
            IsoExonLen = IsoDat[j][4].split(',')
            IsoExonLen.pop()
            IsoExonLen = Covert_Sequence(IsoExonLen)
            
            Dict_Iso_Info[IsoName] = IsoExonNo,IsoExonLab,IsoExonLen
  
        i=i+IsoNo+1
    del Dat


    f_ref = open(FileLen,'r')
    Dict_Exon_Info = {}
    for line in f_ref:
        line = line.rstrip()
        line = line.split('\t')
        gene = line[0]
        ExonLen = line[2].split(',')
        ExonLen.pop()
        ExonLen = Covert_Sequence(ExonLen)
        Dict_Exon_Info[gene] = ExonLen
    f_ref.close()


    f_info = open(ObjectGene,'r')
    Data = [line.rstrip() for line in f_info.readlines()]
    Data = [line.split('\t') for line in Data]
    Data = [line[:2] for line in Data]
    f_info.close()

    parts = 32
    start = 0 
    end = len(Data)
    step = end/parts+1

    job_server = pp.Server()
    job_server.get_ncpus()

    # job_server.set_ncpus(1)
    jobs = []

    for index in xrange(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step,end)
        Datai = Data[starti:endi]
        jobs.append(job_server.submit(Allocate_Location_Workers,(Datai,InputPath,OutputPath,ReadLen,Dict_Gene_Info,Dict_Iso_Info,Dict_Exon_Info,ReadType),globals=globals()))


    for job in jobs:
        job()

def Allocate_Location_Workers(Datai,InputPath,OutputPath,ReadLen,Dict_Gene_Info,Dict_Iso_Info,Dict_Exon_Info,ReadType):

    Len = len(Datai)
    for i in xrange(Len):
        DataLine = Datai[i]

        GeneName = DataLine[0]       
        IsoNo = int(DataLine[1])

        InputFile = os.path.join(InputPath,GeneName)
        OutputFile = os.path.join(OutputPath,GeneName)

        if os.path.isfile(InputFile): 
            Calculate_GeneLocation(GeneName,IsoNo,InputFile,OutputFile,ReadLen,Dict_Gene_Info,Dict_Iso_Info,Dict_Exon_Info,ReadType)
        else:
            pass 



def statisticGeneMap(ObjectGene,OutputPath,FileMap):

    f_ref = open(FileMap,'r')
    Dict_Gene_Map = {}
    for line in f_ref:
        line = line.rstrip()
        line = line.split('\t')
        IsoName = line[1]
        IsoMap = line[2]
        Dict_Gene_Map[IsoName] = IsoMap
    f_ref.close()





    f_info = open(ObjectGene,'r')
    Data = [line.rstrip() for line in f_info.readlines()]
    Data = [line.split('\t') for line in Data]
    f_info.close()


    parts = 32
    start = 0 
    end = len(Data)
    step = end/parts+1

    job_server = pp.Server()
    job_server.get_ncpus()

    # job_server.set_ncpus(4)
    jobs = []

    for index in xrange(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step,end)
        Datai = Data[starti:endi]
        jobs.append(job_server.submit(Allocate_Map_Workers,(Datai,OutputPath,Dict_Gene_Map,),modules=('numpy as np','os','os.path',),globals=globals()))


    for job in jobs:
        job()

def Allocate_Map_Workers(Datai,OutputPath,Dict_Gene_Map):

    Len = len(Datai)
    for i in xrange(Len):
        DataLine = Datai[i]

        GeneName = DataLine[0]       
        IsoNo = int(DataLine[1])

        OutputFile = os.path.join(OutputPath,GeneName)

        IsoNameList = DataLine[3:]

        ModelMap(GeneName,IsoNo,IsoNameList,OutputPath,Dict_Gene_Map)


def ModelMap(GeneName,IsoNo,IsoNameList,OutputPath,Dict_Gene_Map):

    OutputFile = os.path.join(OutputPath,GeneName)

    GeneMap = []
    for i in range(IsoNo):
        IsoName = IsoNameList[i]
        IsoMapList = Dict_Gene_Map[IsoName]
        IsoMapList = IsoMapList.split(',')
        IsoMapList.pop()
        ExonNo = len(IsoMapList)

        IsoMap = []
        for j in range(ExonNo):
            IsoMap.append(IsoMapList[j])
        GeneMap.append(IsoMap)
    
    f_out = open(OutputFile,'w')
    for i in range(ExonNo):
        for j in range(IsoNo):
            f_out.write(str(GeneMap[j][i])+'\t')
        f_out.write('\n')    
    f_out.close()    



def statisticGeneLength(ObjectGene,OutputPath,FileLen):

    f_ref = open(FileLen,'r')
    Dict_Exon_Len = {}
    Dict_Exon_No = {}
    for line in f_ref:    
        line = line.rstrip()
        line = line.split('\t')
        GeneName = line[0].rstrip()
        ExonNo = int(line[1])           
        Dict_Exon_No[GeneName] = ExonNo        
        ExonLen = line[2]
        Dict_Exon_Len[GeneName] = ExonLen
    f_ref.close()


    f_info = open(ObjectGene,'r')
    Data = [line.rstrip() for line in f_info.readlines()]
    Data = [line.split('\t') for line in Data]
    f_info.close()


    parts = 32
    start = 0 
    end = len(Data)
    step = end/parts+1

    job_server = pp.Server()
    job_server.get_ncpus()

    # job_server.set_ncpus(4)
    jobs = []

    for index in xrange(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step,end)
        Datai = Data[starti:endi]
        jobs.append(job_server.submit(Allocate_Length_Workers,(Datai,OutputPath,Dict_Exon_No,Dict_Exon_Len,),modules=('numpy as np','os','os.path',),globals=globals()))


    for job in jobs:
        job()

def Allocate_Length_Workers(Datai,OutputPath,Dict_Exon_No,Dict_Exon_Len):

    Len = len(Datai)
    for i in xrange(Len):
        DataLine = Datai[i]
        GeneName = DataLine[0]  
        OutputFile = os.path.join(OutputPath,GeneName)

        ModelLength(GeneName,OutputPath,Dict_Exon_No,Dict_Exon_Len)


def ModelLength(GeneName,OutputPath,Dict_Exon_No,Dict_Exon_Len):

    ExonLen = Dict_Exon_Len[GeneName]
    ExonLen = ExonLen.split(',')
    ExonNo = Dict_Exon_No[GeneName]

    OutputFile = os.path.join(OutputPath,GeneName)

    f_out = open(OutputFile,'w')
    for i in range(ExonNo):
        f_out.write(ExonLen[i]+'\n')
    f_out.close()


def FunExonLab(Seq):
    Seq = Seq.split(',')
    Seq.pop()
    SeqList = np.array([0 for i in range(len(Seq))])
    for i in range(len(Seq)):
        SeqList[i] = int(Seq[i])
    return SeqList

def FunExonWeight(Seq):
    Seq = FunExonLab(Seq)
    Weight = Seq/sum(Seq+0.0)
    return Weight

      
def statisticGeneData(ObjectGene,InputPath,OutputPath,FileLen):

    f_info = open(ObjectGene,'r')
    Data = [line.rstrip() for line in f_info.readlines()]
    Data = [line.split('\t') for line in Data]
    Data = [line[:3] for line in Data]
    f_info.close()

    f_ref = open(FileLen,'r')
    Dict_Exon_No = {}
    for line in f_ref:    
        line = line.rstrip()
        line = line.split('\t')
        GeneName = line[0].rstrip()
        ExonNo = int(line[1])           
        Dict_Exon_No[GeneName] = ExonNo       
    f_ref.close()


    parts = 32
    start = 0 
    end = len(Data)
    step = end/parts+1

    job_server = pp.Server()
    job_server.get_ncpus()

    # job_server.set_ncpus(1)
    jobs = []

    for index in xrange(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step,end)
        Datai = Data[starti:endi]
        jobs.append(job_server.submit(Allocate_Data_Workers,(Datai,InputPath,OutputPath,Dict_Exon_No),modules=('numpy as np','os','os.path',),globals=globals()))

    for job in jobs:
        job()


def Allocate_Data_Workers(Datai,InputPath,OutputPath,Dict_Exon_No):

    Len = len(Datai)
    for i in xrange(Len):
        DataLine = Datai[i]
        GeneName = DataLine[0]       
        ExonNo = Dict_Exon_No[GeneName]

        ModelData(GeneName,ExonNo,InputPath,OutputPath)


def ModelData(GeneName,ExonNo,InputPath,OutputPath): 

    LaneNo = len(InputPath)
    GeneCount = np.zeros((ExonNo,LaneNo))
    for index in xrange(LaneNo):
    	InputPathi = InputPath[index]
        InputFile = os.path.join(InputPathi, GeneName)

        if os.path.isfile(InputFile):
            f_data = open(InputFile,'r')
            Dict_Read_Info = {}
            for line in f_data:
                line = line.rstrip()
                line = line.split('\t')
                
                Read = line[0]
                if Read not in Dict_Read_Info:
                    Weight = float(line[2])
                    if Weight == 1.0: 
                        ReadLab = FunExonLab(line[5])-1
                        ReadWeight = FunExonWeight(line[6])
                        Dict_Read_Info[Read] = ReadLab,ReadWeight
                    else:
                        pass
                else:
                    pass   
            f_data.close()            
            
            for key in Dict_Read_Info.iterkeys():
                ReadLoc,ReadWeight = Dict_Read_Info[key]
                GeneCount[ReadLoc,index] = GeneCount[ReadLoc,index]+ReadWeight   
        else:
            pass


    OutputFile = os.path.join(OutputPath,GeneName)
    f_out = open(OutputFile,'w')
    for i in range(ExonNo):
        for j in range(LaneNo):
            f_out.write(str(GeneCount[i,j])+'\t')
        f_out.write('\n')    
    f_out.close()
       




