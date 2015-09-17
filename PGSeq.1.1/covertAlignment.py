#! /usr/bin/env python

import pp
import os

def covertAlignment(AnnotationType,AlignmentPaths,AlignmentNames,InputName,ReadType):

	targetDir = os.getcwd()
	referenceDir = os.path.join(targetDir,'Annotation')
	referenceFile = os.path.join(referenceDir,InputName+'.Gene.Label.Isoform')

	f = open(referenceFile,'r')
	DictRef ={}
	for line in f:
		line = line.rstrip()
		line = line.split('\t')
		gene = line[0]
		trans = line[2].rstrip()
		DictRef[trans] = gene
	f.close()


	job_server = pp.Server()
	job_server.get_ncpus()
	# job_server.set_ncpus(4)


	jobs =[]
	covertDataPath = []
	for index in xrange(len(AlignmentPaths)):
		alignPath = AlignmentPaths[index]
		alignName = AlignmentNames[index]

		outputDir = os.path.join(targetDir,'Model.Temp')
		outputDir = os.path.join(outputDir,alignName)
		outputDir = os.path.join(outputDir,'Data.Covert')
		outputFile = os.path.join(outputDir,alignName)

		covertDataPath.append(outputFile)

		jobs.append(job_server.submit(covert,(AnnotationType, alignPath,outputFile,DictRef,ReadType,),globals=globals())) 

	for job in jobs:
		job()


def covert(AnnotationType,InputFile,OutputFile,DictRef,ReadType):
	f = open(InputFile,'r')
	fo = open(OutputFile+'.covert','w')
	if ReadType == 'Single':
		for line in f:
			line = line.rstrip()
			line = line.split('\t')
			read = line[0]
			trans = line[2]
			loc = line[3]


			trans = trans.split('_')
			if AnnotationType in ["Ensembl"]:                   # using annotation from Ensembl
				trans = trans[0]
			elif AnnotationType in ["knownGene",'ensGene']:
				trans = trans[2]
			elif AnnotationType in ["refGene"]:                               # using annotation from UCSC
				trans = trans[2]+'_'+trans[3]

			if trans in DictRef:
				gene = DictRef[trans]
				fo.write(read+'\t'+loc+'\t'+trans+'\t'+gene+'\n')
			else:
				pass

	elif ReadType =='Paired':
		lineNo = 1
		for line in f:
			line = line.rstrip()
			line = line.split('\t')			
			trans = line[2]
			loc = line[3]

			trans = trans.split('_')
			if AnnotationType in ["Ensembl"]:                   # using annotation from Ensembl
				trans = trans[0]
			elif AnnotationType in ["knownGene",'ensGene']:
				trans = trans[2]
			elif AnnotationType in ["refGene"]:                               # using annotation from UCSC
				trans = trans[2]+'_'+trans[3]

			
			if trans in DictRef:
				gene = DictRef[trans]
				if lineNo%2!=0:
					locst = loc
				else:
					read = line[0]
					locnd = loc
					if int(locnd) >= int(locst):
						fo.write(read+'\t'+locst+'\t'+locnd+'\t'+trans+'\t'+gene+'\n')
					else:
						fo.write(read+'\t'+locnd+'\t'+locst+'\t'+trans+'\t'+gene+'\n')

			else:
				pass
			lineNo =lineNo+1


	else:
		pass

	f.close()
	fo.close()

