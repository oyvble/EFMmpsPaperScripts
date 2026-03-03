#rm(list=ls())

#SPECIFY STUDY SETTINGS (used in global script):
AT=10
dynThresh=0

#Include refs
datfold  = "data"
reffile <- paste0(datfold,"/refs.csv") #get file name
refList = sample_tableToList( tableReader( reffile) ) #get sample to analyse
refNames = names(refList)
evidfile <- paste0(datfold,"/evids.csv") #get file name
evidList = sample_tableToList( tableReader( evidfile) ) #get sample to analyse

#Obtain True contributors (listed for each sample)
tab = read.table(file="trueContrTable.csv",sep=";",header=TRUE)
fn = list.files(calcFold,full=T,pattern="DCtables")[1]
sampleNames = names(readRDS(fn))
trueRefList = list()
for(sampleName in sampleNames) {
  #sampleName=sampleNames[17]
  candidates = which(sapply(1:nrow(tab),function(i) grepl(tab[i,1],sampleName)))
  nLetters = nchar(tab[candidates,1])
  maxLetters = which.max(nLetters)
  rowUse = candidates[maxLetters]
  refs = unlist(tab[rowUse,-1])
  refs = refs[refs!=""]
  trueRefList[[sampleName]] = refs
}
sampleNameList = strsplit(sampleNames,"/")
