#THIS SCRIPT CHECKS THE NUMBER OF MISSING MARKERS IN REFERENCE USED IN THE EVALUATION
#rm(list=ls())
maindir = "~\\EFMmpsPaperScripts" 
setwd(maindir)

dsets = list.dirs(recursive = FALSE,full.names = FALSE)#pattern="_")
dsets = dsets[grepl("_",dsets)] #
calcFold = "calcs" #folder where calculations are stored (load from this)

for(ds in dsets) { #for each dataset we run the analysis
#  ds = dsets[1]
  dsdir = paste0(maindir,"\\",ds) #obtain directory of dataset
  setwd(dsdir)
  source("2_createResults.R")
  sampleNames = unlist(sampleNameList)
  sampleNames = intersect(sampleNames,names(evidList))
  #Check number of sites with full ref-registration:
  nMissingLocsRefs = rep(0, length(sampleNames))
  nLocs = rep(0, length(sampleNames))
  for(s in seq_along(sampleNames)) {
#    s=1
    sampleNameVec = sampleNames[s]
    trueRefs = trueRefList[[sampleNameVec[1]]]
    refDatTrue = refList[trueRefs]
    nRefs = length(trueRefs)
    evidListSample = evidList[[s]]
    
    #nAlleles = sapply(evidListSample,function(x) length(x$hdata>=11))
    locs = names(evidListSample)
    nLocs[s] = length(locs)
    nMissingVec = rep(0,nrow=length(locs))
    for(l in seq_along(locs)) {
      loc = locs[l]
      refLocs = lapply(refDatTrue,function(x) unlist(x[[loc]]))
      nMissingVec[l] = any(sapply(refLocs,is.null))
    }
    nMissingLocsRefs[s] = sum(nMissingVec)
  }
  nLocs
  
  range(nMissingLocsRefs)
  range(nLocs)
  #KSNP: "75-88"/"10012-100037"
  #MPS-STR: 0/26
  #MH: 0/100
  #SNP: 2/136
  #CE-STR: 0/15
    
}
  
