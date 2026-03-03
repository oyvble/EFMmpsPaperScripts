#rm(list=ls())
source("0_helpfunctions.R")

#SPECIFY STUDY SETTINGS (used in global script):
AT=11
dynThresh=NULL

#Include refs
allData = getAllData()#onlyREF = TRUE)
refList = allData$refData #Include refs
evidList = allData$evidData

#Obtain True contributors (listed for each sample)
fn = list.files(calcFold,full=T,pattern="DCtables")[1]
sampleNames = names(readRDS(fn))
trueRefList = list()
for(sampleName in sampleNames) trueRefList[[sampleName]] = names(refList)
#compareWithTrueRefs(refList,trueRefList) #go through and create long tables
sampleNameList = strsplit(sampleNames,"/")

