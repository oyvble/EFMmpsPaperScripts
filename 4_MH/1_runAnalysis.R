#rm(list=ls())
#source("0_helpfunctions.R")

methods = c("binned","default","MAE1a","MAE1b","MAE2a","MAE2b")  #methods to consider

#VI option
MAEoptions = list(num_samples = 10000, maxIter = 1000, stepBatch = 50, stopTol = 0.001)

#Model setup
AT = 5  #peak height threshold (static)
lambda=0.0042 #0.01  #dropin hyperparam
modList = list(AT=AT,dynThresh=0.02,imputeDropout = AT/2, pC = 0.05,lambda=0.0042,fst=0.01,
               kit=NULL,DEG=FALSE,BWS=FALSE,FWS=FALSE,SNPmodule=TRUE, nBins=4)

#Read Data
datfold  = "data"
freqfile = list.files(paste0(datfold),pattern="MH_USC_Freq",full=T)
popFreq <- freqImport(freqfile)[[1]]
evidfile <- paste0(datfold,"/evids.csv") #get file name
samplesAll = sample_tableToList( tableReader( evidfile) ) #get sample to analyse
sampleNames <- names(samplesAll) #samplenames

#create a list with samplenames to analyze (also include replicates?)
shortNames <-  sampleNames
shortNamesUnique <-  unique(shortNames) #get sample names (not replicated)
sampleNameList <- NOClist <- list()
for(sind in seq_along(shortNamesUnique)) {
  reps <- sampleNames[shortNamesUnique[sind]==shortNames] #get samples
  for(rep in reps) sampleNameList[[length(sampleNameList)+1]] <- rep
  for(rep in reps) NOClist[[rep]] <- nchar(rep) #depends on sample name
} 




