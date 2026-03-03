#THIS SCRIPT WILL REPRODUCE RESULTS FOUND IN PAPER (THIS IS FOR SUPPLEMENTARY MATERIAL C)
#rm(list=ls())
library(EFMmps)
maindir = "~\\EFMmpsPaperScripts" 
setwd(maindir)

source("1_calculateWithDC_noCond.R")
source("2_compareWithTrueRefs.R") #helpfunction for comparing results
source("3_createStudyResults.R") #helpfunction for creating results (study-wise)
#source("5_createOtherResults.R") #helpfunction for creating results (study-wise)

outfold = "Supplementary" #all results stored here

#Obtain folders with the different datasets
dsets = list.dirs(recursive = FALSE,full.names = FALSE)#pattern="_")
dsets = dsets[grepl("_",dsets)] #

for(ds in dsets) { #for each dataset we run the analysis
#  ds = dsets[3]
  dsdir = paste0(maindir,"\\",ds) #obtain directory of dataset
  setwd(dsdir)
  
  if(ds!="3_MPS-STR") next #only run for this dataset
  source("3_runSupplementary.R") #

  #PERFORM SUPPLEMENTARY CALCULATIONS (stored in Supplementary folder)  
#  methods = "MAE2b_ND" 
  calculateWithDC_noCond(NOClist,sampleNameList,samplesAll,popFreq,modList,MAEoptions,methods,outfold,MAEshrink,MAEshrink2)  #EXECUTING CALCULATIONS

  #CREATE RESULTS (as normal but stored in separate folder)
  if(ds==c("1_KSNP","4_MH")) dynThresh = NULL #deactivate for these two datasets 
  estimateMxWithTrueRefs(sampleNameList,samplesAll,refList,trueRefList,outfold,AT,dynThresh)
  compareDCWithTrueRefs(refList,outfold) #go through and create long tables
  makeMxResults(outfold,outfold) #obtain Mx and PHvar results (fast)
  df = readRDS(file=paste0(outfold,"/ComparedWithTrue_allProbs.RDS"))
  dfSummary = readRDS(file=paste0(outfold,"/ComparedWithTrue_summaryTable.RDS"))
  ruleList =  readRDS(file=paste0(outfold,"/ruleList.RDS"))
  createStudyResults(df, dfSummary, ruleList, outfold) 
  
  #Additional task for STR-MPS:
  if(ds=="3_MPS-STR") {
    source("4_compareDegradParam.R")
    source("5_compareMetrics_againstNoDeg.R")
    source("6_getTimes.R")
  }
}

