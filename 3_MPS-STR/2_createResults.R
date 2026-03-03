  #script to analyze results:
#rm(list=ls())
library(EFMmps)
source("0_helpfunctions.R")

#SPECIFY STUDY SETTINGS (used in global:
AT=30
dynThresh=0

 #Include refs
allData = getAllData()
refList = allData$refData#[1:3] #use only first 3 refs (these are the correct ones)
refNames = names(refList)
evidList = allData$evidData

#Obtain True contributors (listed for each sample)
# tab = read.table(file="trueContrTable.csv",sep=";",header=TRUE)
fn = list.files(calcFold,full=T,pattern="DCtables")[1]
sampleNames = names(readRDS(fn))
trueRefList = list()
for(sampleName in sampleNames) {
  #sampleName=sampleNames[1]
  if(substr(sampleName,1,1)=="2") refs = refNames[1:2] #except for 2-person mixtures
  if(substr(sampleName,1,1)=="3") refs = refNames[1:3] #except for 2-person mixtures
  trueRefList[[sampleName]] = refs
}
sampleNameList = strsplit(sampleNames,"/")


