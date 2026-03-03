#rm(list=ls())
#source("0_helpfunctions.R")

methods = c("binned","default","MAE1a","MAE1b","MAE2a","MAE2b")  #methods to consider

#VI option
MAEoptions = list(num_samples = 10000, maxIter = 1000, stepBatch = 50, stopTol = 0.001)

#Model setup
AT = 10  #peak height threshold
dynThresh=0
modList = list(AT=AT,dynThresh=dynThresh,imputeDropout = AT/2, pC = 0.05,lambda=0.01,fst=0.01,
               kit=NULL,DEG=FALSE,BWS=FALSE,FWS=FALSE,SNPmodule=TRUE, nBins=5)

#Read Data
datfold  = "data"
freqfile = list.files(paste0(datfold,"/Freqs"),pattern="SNP_EUR",full=T)
popFreq <- freqImport(freqfile)[[1]]
evidfile <- paste0(datfold,"/evids.csv") #get file name
samplesAll = sample_tableToList( tableReader( evidfile) ) #get sample to analyse
sampleNames <- names(samplesAll) #samplenames

#create a list with samplenames to analyze (also include replicates?)
shortNames <-  substring(sampleNames,1,nchar(sampleNames)-1) #get sample names (not replicated)
shortNamesUnique <-  unique(shortNames) #get sample names (not replicated)
sampleNameList <- NOClist <- list()
for(sind in seq_along(shortNamesUnique)) {
  reps <- sampleNames[shortNamesUnique[sind]==shortNames] #get samples
  for(rep in reps) sampleNameList[[length(sampleNameList)+1]] <- rep
  for(rep in reps) NOClist[[rep]] <- length(strsplit(rep,"_")[[1]]) #Obtain number of contributor based on name
  sampleNameList[[length(sampleNameList)+1]] <- reps #insert both replicates
} 





