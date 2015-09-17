#!/usr/bin/env python
import os, os.path

def createSubFolder(Path):
	if os.path.exists(Path):
		pass
	else:
		os.mkdir(Path)

def createFolder(Alignment_Folders):

	targetDir = os.getcwd()

	ModelDataPath = os.path.join(targetDir,'PGSeq.Results')
	createSubFolder(ModelDataPath)

	ModelDataPath = os.path.join(targetDir,'Model.Data')
	createSubFolder(ModelDataPath)

	createSubFolder(os.path.join(ModelDataPath,'Gene.Map'))
	createSubFolder(os.path.join(ModelDataPath,'Gene.Data'))
	createSubFolder(os.path.join(ModelDataPath,'Gene.Length'))

	ModelTempPath = os.path.join(targetDir,'Model.Temp')
	createSubFolder(ModelTempPath)

	for i in xrange(len(Alignment_Folders)):
		AlignFolder = Alignment_Folders[i]
		AlignFolderPath = os.path.join(ModelTempPath,AlignFolder)

		createSubFolder(AlignFolderPath)

		createSubFolder(os.path.join(AlignFolderPath,'Gene.LocAbs'))
		createSubFolder(os.path.join(AlignFolderPath,'Gene.LocRel'))
		createSubFolder(os.path.join(AlignFolderPath,'Data.Prob'))
		createSubFolder(os.path.join(AlignFolderPath,'Data.Covert'))


def deleteTempFolder():	
    '''delete files and folders''' 
    targetDir = os.getcwd()
    Path = os.path.join(targetDir,'Model.Temp')
    deleteSubFolder(Path)

def deleteSubFolder(Path):
    if os.path.isfile(Path):  
        try:  
            os.remove(Path)  
        except:  
            pass 
    elif os.path.isdir(Path):  
        for item in os.listdir(Path):  
            itemPath=os.path.join(Path,item)  
            deleteSubFolder(itemPath)  
        try:  
            os.rmdir(Path)  
        except:  
            pass 


def deleteDataFolder():	
    '''delete files and folders''' 
    targetDir = os.getcwd()
    Path = os.path.join(targetDir,'Model.Data')
    deleteSubFolder(Path)













