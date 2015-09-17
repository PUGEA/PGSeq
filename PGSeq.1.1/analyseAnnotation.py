# /usr/bin/env python

import os, os.path
import sys


def Create_Dic_IsoNum(table):
    Dic_IsoNum = {}
    for i in range(len(table)):
        GeneName = table[i][0]        
        if GeneName not in Dic_IsoNum:
           Dic_IsoNum[GeneName] = 1
        else:
           Dic_IsoNum[GeneName] = Dic_IsoNum[GeneName] + 1
    return  Dic_IsoNum


def SelectedGenesInfo(targetGene,Prefix,InputDir):

	allGeneInfo = os.path.join(InputDir,Prefix)
	allGeneInfoFile = allGeneInfo+'.Gene.Info'  
	f=open(allGeneInfoFile,'r')
	DictGeneInfo = {}
	for line in f:
		line = line.rstrip()
		line = line.split('\t')
		gene = line[0]
		info = line[1:]
		DictGeneInfo[gene] = info
	f.close()

	OutputFile = os.path.join(InputDir,'Selected.Genes.Info')

	fsel = open(targetGene,'r')
	fout = open(OutputFile,'w')
	for line in fsel:
		line = line.rstrip()
		line = line.split('\t')
		gene = line[0]

		if gene in DictGeneInfo:
			info = DictGeneInfo[gene]
			fout.write(gene+'\t')
			for i in xrange(len(info)):
				fout.write(info[i]+'\t')
			fout.write('\n')
		else:
			print gene+' is not in the annnotation files.'

	fsel.close()
	fout.close()

	return OutputFile


def CovertAnnotation(AnnotationType,AnnotationFile,OutputFile,OutputDir):

	if AnnotationType in ["ensGene","refGene"]:
		f_ens = open(AnnotationFile,'r')
		covertFile = os.path.join(OutputDir,OutputFile+'.covertFormat')
		f_flat = open(covertFile,'w')
		f_ens.readline()
		for line in f_ens:
		    line = line.rstrip()
		    line = line.split('\t')

		    trans = line[1]
		    chrom,strand = line[2:4]
		    txStart,txEnd,cdsStart,cdsEnd = line[4:8]
		    exonCount,exonStarts,exonEnds = line[8:11]
		    gene = line[12]

		    f_flat.write(gene+'\t'+trans+'\t'+chrom+'\t'+strand+'\t'+txStart+'\t'+txEnd+'\t'+cdsStart+'\t'+cdsEnd+'\t'+exonCount+'\t')
		    exonStarts = exonStarts.split(',')
		    for i in range(int(exonCount)):
		        f_flat.write(exonStarts[i]+',')
		    f_flat.write('\t')

		    exonEnds = exonEnds.split(',')
		    for i in range(int(exonCount)):
		        f_flat.write(exonEnds[i]+',')
		    f_flat.write('\n')

		f_ens.close()
		f_flat.close()

	elif AnnotationType in ["knownGene"]:

		AnnotationFileList = AnnotationFile.split(',')
		if len(AnnotationFileList)==2:
			pass
		else:
			print 'This programming needs to knownGene.txt and knwonIsoforms.txt for "knownGene" option.'
			sys.exit()


		f = open(AnnotationFileList[0],'r')
		f.readline()
		isoInfo ={}
		for line in f:
			line = line.rstrip()
			line = line.split('\t')
			isof = line[0].strip()
			info = line[1:]
			isoInfo[isof]=info
		f.close()


		f = open(AnnotationFileList[1],'r')
		covertFile = os.path.join(OutputDir,OutputFile+'.covertFormat')
		f_flat = open(covertFile,'w')
		f.readline()
		for line in f:
			line = line.rstrip()
			line = line.split('\t')
			group, isof = line
			isof = isof.strip()
			group = group.strip()
			info = isoInfo[isof]

			exno = int(info[6])
			star = info[7].split(',')
			end = info[8].split(',')
			refgene = info[9]

			f_flat.write(group+'\t'+isof+'\t')
			for i in range(7):
				f_flat.write(info[i]+'\t')			
			for i in range(exno):
				f_flat.write(star[i]+',')
			f_flat.write('\t')
			for i in range(exno):
				f_flat.write(end[i]+',')
			f_flat.write('\n')
		f.close()
		f_flat.close()


	elif AnnotationType in ["Ensembl"]:
		Dict_seqTran = {}

		f= open(AnnotationFile,'r')
		for line in f:
			line = line.rstrip()
			line = line.split('\t')

			info = line[8]
			info = info.split('"')
			trans = info[3]

			if trans not in Dict_seqTran:
				Dict_seqTran[trans] = trans
			else:
				pass



		f = open(AnnotationFile,'r')
		dict_gen = {}
		dict_chr = {}
		dict_tra_start ={}
		dict_tra_end = {}
		dict_tra_strand ={}
		for line in f:
			line = line.rstrip()
			line = line.split('\t')

			chrN = line[0]
			exon = line[2]
			start = int(line[3])-1
			start = str(start)
			end = line[4]
			strand = line[6]

			info = line[8]
			info = info.split('"')
			gene = info[1]
			trans = info[3]
			# exon_lab = info[5]
			# gene_ucsc = info[7]
			# trans_ucsc = info[11]

			if exon == "exon":
				if gene not in dict_gen:
					dict_gen[gene] = [trans]
				else:
					temp = dict_gen[gene]
					if trans in temp:
						pass 
					else:
						temp.append(trans)
						dict_gen[gene] = temp

				if trans not in dict_tra_start:
					dict_chr[trans] = chrN
					dict_tra_strand[trans] = strand
					dict_tra_start[trans] = [start]
					dict_tra_end[trans] = [end]
				else:					
				    tempstart = dict_tra_start[trans]
				    tempstart.append(start)
				    dict_tra_start[trans] = tempstart

				    tempend = dict_tra_end[trans]
				    tempend.append(end)
				    dict_tra_end[trans] = tempend
			else:
				pass
		f.close()

		covertFile = os.path.join(OutputDir,OutputFile+'.covertFormat')
		f_flat = open(covertFile,'w')
		for gene in dict_gen.iterkeys():
			translist = dict_gen[gene]
			for i in range(len(translist)):
				trans = translist[i]
				if trans in Dict_seqTran:
					chrN = dict_chr[trans]
					strand = dict_tra_strand[trans]
					start = dict_tra_start[trans]
					end = dict_tra_end[trans]

					f_flat.write(gene+'\t'+trans+'\t'+chrN+'\t'+strand+'\t')
					for j in range(4):
						f_flat.write('*'+'\t')
					f_flat.write(str(len(start))+'\t')



					if strand == '+':
						for j in range(len(start)):
							f_flat.write(start[j]+',')
						f_flat.write('\t')
						for j in range(len(end)):
							f_flat.write(end[j]+',')
						f_flat.write('\n')
					else:
						start.reverse()
						end.reverse()
						for j in range(len(start)):
							f_flat.write(start[j]+',')
						f_flat.write('\t')
						for j in range(len(end)):
							f_flat.write(end[j]+',')
						f_flat.write('\n')    
				else:
					pass

		f_flat.close()
	else:
		pass


def AnalyseAnnotation(AnnotationFile,OutputFile,OutputDir):


	# print 'Loading the covert  table.'
	covertFile = os.path.join(OutputDir,OutputFile+'.covertFormat')
	f_in = open(covertFile,'r')
	Dat = [line.strip() for line in f_in.readlines()]
	Dat = [line.split('\t') for line in Dat]
	Dat.sort()
	f_in.close()


	##------------------------------------------------------------##
	## Check duplicates transcripts, Skipped.
	##------------------------------------------------------------##
	# print '      Checking duplicate transcripts.'

	Dic_IsoNum = {}
	Dic_IsoNum = Create_Dic_IsoNum(Dat)
	CountDupl = 0 
	i = 0
	NoDuplDat = []           # save the new data without duplicate transcripts 
	while i < len(Dat):
	    DatTemp = Dat[i]
	    GeneName = DatTemp[0]
	    IsoName = DatTemp[1]
	    IsoNum = Dic_IsoNum[GeneName]


	      
	    if IsoNum == 1:
	       NoDuplDat.append(DatTemp)
	       i=i+1
	    else:
	       NoDuplDat.append(DatTemp)
	       IsoNameList = [IsoName]
	       DatTempx = Dat[i+1:i+IsoNum]
	       for j in range(len(DatTempx)):
	           Temp = DatTempx[j]
	           IsoG = DatTempx[j][1]        
	           if IsoG not in IsoNameList:
	              IsoNameList.append(IsoG)
	              NoDuplDat.append(Temp)              
	           else: 
	              CountDupl =CountDupl+1
	#              print '%s \t  %s ' %(GeneName, IsoG)
	              if(CountDupl<=10):
		              print '      %s have duplcate %s transctipts,%s chromosome skipped' %(GeneName, IsoG,DatTempx[j][2])
	              #pass
	       i = i+IsoNum 
	if CountDupl == 0:
		pass
	else:
		print '      %d duplicate transcripts, skipped' %(CountDupl)

	del Dat,Dic_IsoNum, DatTemp,DatTempx
	##------------------------------------------------------------##
	## Check transcript consistency, Chr or Strand Removed 
	##------------------------------------------------------------##
	# print '      Checking transcript consistency.'
	Dic_NoDuplIsoNum = {}
	Dic_NoDuplIsoNum = Create_Dic_IsoNum(NoDuplDat)


	i = 0
	ConsistDat = []
	while i < len(NoDuplDat):
	    GeneName = NoDuplDat[i][0]
	    IsoNum = Dic_NoDuplIsoNum[GeneName]

	    if IsoNum == 1:
	       ConsistDat.append(NoDuplDat[i])
	       i = i + 1
	    else:
	       ConsistDat.append(NoDuplDat[i])
	       Chr = NoDuplDat[i][2]
	       Strand = NoDuplDat[i][3]
	       ChrList = [Chr]
	       StrandList = [Strand]

	       DatTemp = NoDuplDat[i+1:i+IsoNum]
	       for j in range(len(DatTemp)):
	           Temp = DatTemp[j]
	           IsoName = DatTemp[j][1] 
	           ChrG = DatTemp[j][2]
	           StrandG = DatTemp[j][3]
	           ChrList.append(ChrG)
	           StrandList.append(StrandG)
	           # ConsistDat.append(Temp)
	           if ChrG in ChrList and StrandG in StrandList:
	              ChrList.append(ChrG)
	              StrandList.append(StrandG)
	              ConsistDat.append(Temp)
	           else:
	              if ChrG not in ChrList:
	                  print '      %s %s Chr is not consistency, removed' %(GeneName,IsoName)
	              else:
	                  print '      %s %s Strand is not consistency, removed' %(GeneName,IsoName)


	       i = i+IsoNum           
	del NoDuplDat,Dic_NoDuplIsoNum

	##------------------------------------------------------------##
	## Check transcript coordinates, Exon and intron union 
	##------------------------------------------------------------##
	# print '      Checking transcript coordinates.'
	i = 0
	FinalDat = []
	BadExon =0
	BadIntron=0

	while i < len(ConsistDat):
	    GeneName = ConsistDat[i][0]
	    IsoName = ConsistDat[i][1]
	    ExonNum = int(ConsistDat[i][8])
	    StartTemp = ConsistDat[i][9].split(',')
	    EndTemp = ConsistDat[i][10].split(',')

	    ExonStart =[]
	    ExonEnd = []
	    for j in range(ExonNum):
	        ExonStart.append(int(StartTemp[j]))
	        ExonEnd.append(int(EndTemp[j]))
	    ExonStartL = ExonStart[:]
	    ExonEndL = ExonEnd[:]    

	    GoodExon = 1
	    ExonNumFinal = ExonNum
	    for j in range(ExonNum):
	        if(not (ExonStart[j]>=0 and ExonEnd[j] >= ExonStart[j] and (j==0 or ExonStart[j] >= ExonEnd[j-1]))):
	             GoodExon = 0
	             break
	        if(ExonStart[j] == ExonEnd[j]):
	             ExonStartL.remove(ExonStart[j])
	             ExonEndL.remove(ExonEnd[j])
	             ExonNumFinal = ExonNumFinal - 1

	             # BadExon = BadExon+1
	             # if BadExon <=10:
	             # 	print 'Warning: %s %s exon length is zero, removed' %(GeneName,IsoName)
	             # else:
	             # 	pass

	             continue

	        if(j>0 and ExonEnd[j-1] == ExonStart[j]):
	             ExonStartL.remove(ExonStart[j])
	             ExonEndL.remove(ExonEnd[j-1])
	             ExonNumFinal = ExonNumFinal -1

	            #  BadIntron = BadIntron +1
	            #  if BadIntron<=10:
	            #  	print 'Warning: %s %s intron length is zero, merged' %(GeneName,IsoName)
	            #  else:
 		         	# pass
 		         	
	             continue

	    if(GoodExon == 0):
	        print 'Warning: %s %s is incoordinates, removed' %(GeneName,IsoName)
	    else:
	        DatTemp = []
	        DatTemp = ConsistDat[i][:8]
	        DatTemp.append(ExonNumFinal)
	        DatTemp.append(ExonStartL)
	        DatTemp.append(ExonEndL)

	        FinalDat.append(DatTemp) 

	    i = i + 1


	##------------------------------------------------------------##
	## subexon split
	##------------------------------------------------------------##
	# print '      Spliting Exon without overlap.'

	Dic_FinalIsoNum = {}
	Dic_FinalIsoNum = Create_Dic_IsoNum(ConsistDat)

	SubExonDat = []
	i = 0
	while i < len(FinalDat):
	     DatLine = FinalDat[i]
	     GeneName = DatLine[0]
	     IsoNum = Dic_FinalIsoNum[GeneName]
	 
	     if IsoNum == 1:
	        IsoName = DatLine[1]
	        ExonNum = DatLine[8]
	        ExonIndex = [1]*int(ExonNum)
	         
	        GeneTemp = [GeneName]+DatLine[2:4]+[IsoNum]+DatLine[8:11]
	        SubExonDat.append(GeneTemp)
	        IsoTemp = [IsoName]+DatLine[2:4]+[1]+[ExonNum]+[ExonIndex]
	        SubExonDat.append(IsoTemp)

	        i=i+1
	     else:
	        ExonLocS = []
	        ExonLocE = []
	        ExonExist = []        
	        DatMat = FinalDat[i:i+IsoNum]
	        for j in range(IsoNum):
	            ExonLocS = ExonLocS + DatMat[j][9]
	            ExonLocE = ExonLocE + DatMat[j][10]


	        for j in range(len(ExonLocS)):
	            ExonExist.append([ExonLocS[j],ExonLocE[j]])
	        ExonLoc = ExonLocS + ExonLocE

	        ExonUniqLoc = []
	        for item in set(ExonLoc):
	            ExonUniqLoc.append(item)
	        ExonUniqLoc.sort()     

	        ExonCell = []
	        for j in range(len(ExonUniqLoc)-1):
	            ExonCell.append([ExonUniqLoc[j],ExonUniqLoc[j+1]])

	        ExonSplit = []
	        for j in range(len(ExonCell)):
	            LocSC,LocEC = ExonCell[j]
	            for k in range(len(ExonExist)):
	                LocSEx, LocEEx = ExonExist[k]
	                if(LocSC >= LocSEx and LocEC <= LocEEx):
	                      ExonSplit.append([LocSC,LocEC])
	                      break
	                else: 
	                      pass 

	        LocSTemp = []
	        LocETemp = []
	        for j in range(len(ExonSplit)):
	            LocSTemp.append(ExonSplit[j][0])
	            LocETemp.append(ExonSplit[j][1])
	        GeneTemp = [GeneName]+DatLine[2:4]+[IsoNum]+[len(ExonSplit)]+[LocSTemp]+[LocETemp]    
	        SubExonDat.append(GeneTemp)      



	        for j in range(IsoNum):
	            IsoName = DatMat[j][1]
	            ExonNum = DatMat[j][8]
	            LocStart = DatMat[j][9]
	            LocEnd = DatMat[j][10] 
	            LocIndex = []

	            for k in range(len(ExonSplit)):
	                LocSSp = ExonSplit[k][0]
	                LocESp = ExonSplit[k][1]
	                Flag = 0 
	                for l in range(ExonNum):
	                    if(LocSSp>=LocStart[l] and LocESp<=LocEnd[l]):
	                       Flag = 1
	                    else:
	                       pass
	                LocIndex.append(Flag)

	            IsoTemp = [IsoName]+DatMat[j][2:4]+[j+1]+[ExonNum]+[LocIndex]
	            SubExonDat.append(IsoTemp)
	  
	   
	        i=i+IsoNum





	##------------------------------------------------------------##
	## Write the subexon file 
	#------------------------------------------------------------##
	# print 'Writing the ExonSplit file.'

	ExonSplitFile = os.path.join(OutputDir,OutputFile+'.ExonSplit')

	f_out = open(ExonSplitFile,'w')

	i = 0
	while i < len(SubExonDat):
	    GeneDat = SubExonDat[i]
	    GeneName = GeneDat[0]
	    IsoNum = GeneDat[3]
	         
	    StartTemp = GeneDat[5]
	    EndTemp = GeneDat[6]   
	    ExonStart = []
	    ExonEnd = []
	    for j in range(len(StartTemp)):
	        ExonStart.append(str(StartTemp[j]))
	        ExonEnd.append(str(EndTemp[j]))    
	    ExonStart = ','.join(ExonStart)
	    ExonEnd = ','.join(ExonEnd)  

	    for j in range(len(GeneDat)-2):
	        f_out.writelines(str(GeneDat[j])+'\t')  
	    f_out.writelines(ExonStart+'\t')
	    f_out.writelines(ExonEnd+'\n')            

	    IsoDat = SubExonDat[i+1:i+IsoNum+1]

	    for j in range(IsoNum):
	        IsoDatT = IsoDat[j]
	        for k in range(len(IsoDatT)-1):
	            f_out.writelines(str(IsoDatT[k])+'\t')             
	        IndexT = IsoDatT[5]
	        IndexList = []
	        for k in range(len(IndexT)):
	            IndexList.append(str(IndexT[k]))
	        IndexList = ','.join(IndexList)
	        f_out.writelines(IndexList+'\n')
	    i = i+IsoNum+1 

	f_out.close()

	#-------------------------------------------------------------##          


	f_in= open(ExonSplitFile,'r')
	Dat = [line.strip() for line in f_in.readlines()]
	Dat = [line.split('\t') for line in Dat]
	f_in.close()

	f_t1 = open(os.path.join(OutputDir,OutputFile+'.Exon.Split'),'w')
	f_t2 = open(os.path.join(OutputDir,OutputFile+'.Exon.Len'),'w')
	f_t3 = open(os.path.join(OutputDir,OutputFile+'.Gene.Label.Isoform'),'w')
	f_t4 = open(os.path.join(OutputDir,OutputFile+'.Gene.Map.Isoform'),'w')
	f_t5 = open(os.path.join(OutputDir,OutputFile+'.Gene.Len'),'w')
	f_t6 = open(os.path.join(OutputDir,OutputFile+'.Gene.Info'),'w')

	LineNo = 0
	# Select and process a gene
	while LineNo < len(Dat):
	      GeneData = Dat[LineNo]
	      GeneName = GeneData[0]
	      Dir = GeneData[2]
	      IsoNo = int(GeneData[3])
	      ExoNo = int(GeneData[4])
	      ExonStartLocList = GeneData[5].split(',')
	      ExonEndLocList = GeneData[6].split(',')
	      
	      if Dir == '-':
	            ExonStartLocList.reverse()
	            ExonEndLocList.reverse()
	      
	      GeneExonLabList = []
	      GeneExonLenList = []
	      GeneExonIterLenList = []
	      GeneExonLenSum = 0
	      for i in range(ExoNo):  
	            ExonStart = int(ExonStartLocList[i])
	            ExonEnd = int(ExonEndLocList[i])
	            ExonLen = ExonEnd - ExonStart
	            GeneExonLenSum = GeneExonLenSum + ExonLen
	            GeneExonLabList.append(i+1)
	            GeneExonLenList.append(ExonLen)
	            GeneExonIterLenList.append(GeneExonLenSum)
	            
	 #----------  output  start     -------------------------------  
	      
	      f_t1.write(GeneName+'\t'+str(IsoNo)+'\t'+str(ExoNo)+'\t')
	      for i in range(ExoNo):
	            f_t1.write(str(GeneExonLabList[i])+',')
	      f_t1.write('\t')
	      for i in range(ExoNo):
	            f_t1.write(str(GeneExonIterLenList[i])+',')
	      f_t1.write('\n')
	      
	      f_t2.write(GeneName+'\t'+str(ExoNo)+'\t')
	      for i in range(ExoNo):
	            f_t2.write(str(GeneExonLenList[i])+',')
	      f_t2.write('\n')

	      f_t5.write(GeneName+'\t'+str(IsoNo)+'\t'+str(ExoNo)+'\t')
	      for i in range(ExoNo):
	            f_t5.write(str(GeneExonIterLenList[i])+',')
	      f_t5.write('\n')


	      f_t6.write(GeneName+'\t'+str(IsoNo)+'\t'+str(GeneExonIterLenList[-1])+'\t')



	#----------  output   end    -------------------------------   

	      IsoDataList = Dat[LineNo+1:LineNo+IsoNo+1]
	      for i in range(IsoNo):
	            IsoName = IsoDataList[i][0]
	            IsoLab = IsoDataList[i][3]
	            IsoSplit = IsoDataList[i][5].split(',') 
	            
	            if Dir =='-':
	                  IsoSplit.reverse()
	            
	            
	            IsoExonLabList = []
	            IsoExonLenList = []
	            IsoExonIterLenList = []
	            IsoExonLenSum = 0
	            for j in range(ExoNo):
	                  if IsoSplit[j] == '0':
	                        pass
	                  else:
	                        IsoExonLabList.append(GeneExonLabList[j])
	                        IsoExonLenList.append(GeneExonLenList[j])
	                        IsoExonLenSum = IsoExonLenSum + GeneExonLenList[j]
	                        IsoExonIterLenList.append(IsoExonLenSum)
	                        
	            IsoExoNo = len(IsoExonLabList)
	            f_t1.write(GeneName+'\t'+IsoName+'\t'+str(IsoExoNo)+'\t')
	            for j in range(IsoExoNo):
	                  f_t1.write(str(IsoExonLabList[j])+',')
	            f_t1.write('\t')
	            for j in range(IsoExoNo):
	                  f_t1.write(str(IsoExonIterLenList[j])+',')
	            f_t1.write('\n')
	                        

	            f_t3.write(GeneName+'\t'+IsoLab+'\t'+IsoName+'\n')

	            f_t4.write(GeneName+'\t'+IsoName+'\t')
	            for j in range(ExoNo):
	                  f_t4.write(IsoSplit[j]+',')
	            f_t4.write('\n')

	            f_t6.write(IsoName+'\t')

	      f_t6.write('\n')
	      LineNo = LineNo+IsoNo+1

	f_t1.close()
	f_t2.close()
	f_t3.close()
	f_t4.close()
	f_t5.close()
	f_t6.close()

	os.remove(covertFile)
	os.remove(ExonSplitFile)

