  #script to analyze results:
#rm(list=ls())
library(EFMmps)
source("0_helpfunctions.R")

#SPECIFY STUDY SETTINGS (used in global:
AT=50
dynThresh=0

 #Include refs
allData = getAllData()
refList = allData$refData#[1:3] #use only first 3 refs (these are the correct ones)
refNames = names(refList)
evidList = allData$evidData

#Obtain True contributors (listed for each sample)
# tab = read.table(file="trueContrTable.csv",sep=";",header=TRUE)
fn = list.files(calcFold,full=T,pattern="DCtables")[1] #read one of the result files
sampleNames = names(readRDS(fn))
trueRefList = list()
for(sampleName in names(evidList)) {
#  sampleName=sampleNames[1]
  sampleNameSplitted = strsplit(sampleName,"\\.")[[1]]
  contrID = sampleNameSplitted[1]
  refGroup = substr(refNames,1,nchar(refNames)-1)
  refs = refNames[which(refGroup==contrID)] #Get the correct refs (3p by default)
  
  shortName = paste0(sampleNameSplitted[1:2],collapse=".") #obtain the short name (not replicated)
  if(shortName%in%c("0.5","0.9","0.24","0.28")) refs = refs[1:2] #these are 2p
  trueRefList[[sampleName]] = refs
}

sampleNames2 = sampleNames #make a copyt
isCombinedInds = grep("/",sampleNames)
for(s in isCombinedInds) {
  #    s=isCombinedInds[1]
  sampleList = strsplit(sampleNames[s],"/")[[1]]
  sampleNames2[sampleNames%in%sampleList] = NA
}
sampleNames = na.omit(sampleNames2)
sampleNameList = strsplit(sampleNames,"/")


