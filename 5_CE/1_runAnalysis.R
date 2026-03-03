#rm(list=ls())
source("0_helpfunctions.R")

methods = c("default","MAE1a","MAE1b","MAE2a","MAE2b")  #methods to consider

#VI option
MAEoptions = list(num_samples = 10000, maxIter = 1000, stepBatch = 50, stopTol = 0.001)

#Model setup
AT = 30  #peak height threshold
modList = list(AT=AT,dynThresh=0,imputeDropout = AT/2, pC = 0.00036,lambda=0.02,fst=0.01,
               kit="NGM",DEG=TRUE,BWS=FALSE,FWS=FALSE,SNPmodule=FALSE)

#Read Data
allData = getAllData()
popFreq <- allData$popFreq
samplesAll = allData$evidData
sampleNameListAll = allData$sampleNameList #get all samples (also combined)

#Obtain number of NOC to assume and samples to consider
isCombinedSamples = which(sapply(sampleNameListAll,function(x) length(x)>1))
allCombinedSamples = unlist(sampleNameListAll[isCombinedSamples])
NOClist <- sampleNameList <- list()
for(sind in seq_along(sampleNameListAll)) { #for each samples in list
  #  sind = 1
  NOC = 3
  sampleNameVec = sampleNameListAll[[sind]]

  sampleNameReps <- paste0(sampleNameVec,collapse="/")
  shortName = paste0(strsplit(sampleNameReps,"\\.")[[1]][1:2],collapse=".")
  if(shortName%in%c("0.5","0.9","0.24","0.28")) NOC=2 #these are 2p mixtures
  NOClist[[sampleNameVec[1]]] = NOC
  
  if(length(sampleNameVec)==1 && sampleNameVec%in%allCombinedSamples) next #skip if sample given as combined
  sampleNameList[[length(sampleNameList)+1]] = sampleNameVec
}
  



