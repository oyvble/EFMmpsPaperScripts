#THIS SCRIPT WILL REPRODUCE RESULTS FOUND IN PAPER
#rm(list=ls())
library(EFMmps)
maindir = "~\\EFMmpsPaperScripts" 
setwd(maindir)
source("2_compareWithTrueRefs.R") #helpfunction for comparing results
source("3_createStudyResults.R") #helpfunction for creating results (study-wise)

#Obtain folders with the different datasets
dsets = list.dirs(recursive = FALSE,full.names = FALSE)#pattern="_")
dsets = dsets[grepl("_",dsets)] #

for(ds in rev(dsets)) { #for each dataset we run the analysis
#  ds = dsets[3]
  dsdir = paste0(maindir,"\\",ds) #obtain directory of dataset
  setwd(dsdir)
  print(ds)
  
  #Perform comparison
  calcFold = "calcs" #folder where calculations are stored (load from this)
  outFold = "Results" #folder where results are created (dump here)
  dir.create(outFold)
  #if(0) {
  source("2_createResults.R") #obtain necessary variables from study

  #Following steps are excecuted after defining necessary variables in study:
  
  #Step 1) Estimate Mx with true refs assumed
  #estimateMxWithTrueRefs: Estimating Mx conditioning on true references
  estimateMxWithTrueRefs(sampleNameList,evidList,refList,trueRefList,calcFold,AT,dynThresh)
  
  #Step 2) Compare DC results with true refs (with Mx details)
  #compareDCWithTrueRefs: Comparing DC results with true refs (uses Mx details in prev. step to map against unknown)
  compareDCWithTrueRefs(refList,calcFold) #go through and create long tables
  
  #Step 3) Create Mx and PHvar results (per sample)
  #makeMxResults: obtain Mx and PHvar results per sample (detailed results)
  makeMxResults(calcFold,outFold) #obtain Mx and PHvar results (fast)
  #}
  
  #Step 4) PRODUCING RESULTS FOR PAPER based on stored objects
  df = readRDS(file=paste0(calcFold,"/ComparedWithTrue_allProbs.RDS"))
  dfSummary = readRDS(file=paste0(calcFold,"/ComparedWithTrue_summaryTable.RDS"))
  ruleList =  readRDS(file=paste0(calcFold,"/ruleList.RDS"))
  
  #Create tables and figures for study
  createStudyResults(df, dfSummary, ruleList, outFold) 
}
