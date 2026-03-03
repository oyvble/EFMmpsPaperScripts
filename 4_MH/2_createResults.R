#rm(list=ls())

#SPECIFY STUDY SETTINGS (used in global script):
AT=5
dynThresh=0

#Include refs
datfold  = "data"
reffile <- paste0(datfold,"/refs.csv") #get file name
refList = sample_tableToList( tableReader( reffile) ) #get sample to analyse
refNames = names(refList)
evidfile <- paste0(datfold,"/evids.csv") #get file name
evidList = sample_tableToList( tableReader( evidfile) ) #get sample to analyse

#Obtain True contributors (listed for each sample)
#  tab = read.table(file="trueContrTable.csv",sep=";",header=TRUE)
fn = list.files(calcFold,full=T,pattern="DCtables")[1] #get sample names
sampleNames = names(readRDS(fn))
trueRefList = list()
for(sampleName in sampleNames) {
  #sampleName=sampleNames[1]
  refs = refNames #all 3 are true donors
  if(nchar(sampleName)<3) refs = refNames[1:2] #except for 2-person mixtures
  trueRefList[[sampleName]] = refs
}
sampleNameList = strsplit(sampleNames,"/")


