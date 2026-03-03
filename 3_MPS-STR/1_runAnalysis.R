#rm(list=ls())
source("0_helpfunctions.R")

methods = c("default","MAE1a","MAE1b","MAE2a","MAE2b")  #methods to consider

#VI option
MAEoptions = list(num_samples = 10000, maxIter = 1000, stepBatch = 50, stopTol = 0.001)

#Model setup
AT = 30  #peak height threshold
modList = list(AT=AT,dynThresh=0,imputeDropout = AT/2, pC = 0.05,lambda=0.01,fst=0.01,
               kit="ForenSeq",DEG=TRUE,BWS=TRUE,FWS=FALSE,SNPmodule=FALSE)

#Read Data
allData = getAllData(filtertype="static_withStutters",format="LUS+", AT=AT, skipEmpty=FALSE,datfold = "data")
popFreq <- allData$popFreq
samplesAll = allData$evidData
sampleNames <- names(samplesAll) #samplenames
sampleNames = sampleNames[substr(sampleNames,1,1)!="4"] #not 4 contr

#create a list with samplenames to analyze (also include replicates?)
shortNames <-  substring(sampleNames,3,nchar(sampleNames)) #get sample names (not replicated)
shortNamesUnique <-  unique(shortNames) #get sample names (not replicated)
sampleNameList <- NOClist <- list()
for(sind in seq_along(shortNamesUnique)) {
  reps <- sampleNames[shortNamesUnique[sind]==shortNames] #get samples
  for(rep in reps) sampleNameList[[length(sampleNameList)+1]] <- rep
  for(rep in reps) NOClist[[rep]] <- as.integer(substr(rep,1,1)) #Obtain number of contributor based on name
  sampleNameList[[length(sampleNameList)+1]] <- reps #insert both replicates
} 





