#THIS SCRIPT WILL REPRODUCE RESULTS FOUND IN PAPER
library(EFMmps)
maindir = "~\\EFMmpsPaperScripts" 
setwd(maindir)
source("1_calculateWithDC_noCond.R")

#Obtain folders with the different datasets
dsets = list.dirs(recursive = FALSE,full.names = FALSE)#pattern="_")
dsets = dsets[grepl("_",dsets)] #

for(ds in rev(dsets)) { #for each dataset we run the analysis
#  ds = dsets[5]
  dsdir = paste0(maindir,"\\",ds) #obtain directory of dataset
  setwd(dsdir)
  source("1_runAnalysis.R")
 
  #methods = "MAE2a"  #force this
  #EXECUTING CALCULATIONS (variabels are defined globally)
  #methods = c("MAE1a","MAE1b")#,"MAE2a","MAE2b")  #methods to consider
  calculateWithDC_noCond(NOClist,sampleNameList,samplesAll,popFreq,modList,MAEoptions,methods)  
}
