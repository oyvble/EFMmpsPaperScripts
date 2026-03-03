
#type = "evid"

getAllData = function(datfold="data") {
  locrm = "AMEL"
  databaseFile = paste0(datfold,"\\NGM_Holland2.csv") #The allele frequency file
  popFreq <- euroformix::freqImport(databaseFile)[[1]] #import population freqs.

  reffiles <- list.files(paste0(datfold,"\\References"),recursive = FALSE,full=T) #obtain all folders (each single is a sample)
  #refNames = sapply(strsplit(basename(reffiles),"\\."),function(x) x[1])
  
  refData = list()
  for(r in seq_along(reffiles)) {
    #  r=1
    refObj = sample_tableToList(tableReader(reffiles[r])) #get evid element
    refData[names(refObj)] = refObj
  }
  length(refData)
  
  #Get files
  evidfiles <- list.files(paste0(datfold,"\\SamplesStutterFilter"),recursive = FALSE,full=T) #obtain all folders (each single is a sample)
  evidNames = sapply(strsplit(basename(evidfiles),"_"),function(x) x[1])
  

  evidData = NULL
  sampleNameList = list()
  for(e in seq_along(evidfiles)) {
#    e=1
    evidfile = evidfiles[e]
    evidTab = tableReader(evidfile)
    evidTab = evidTab[!toupper(evidTab$Marker)%in%locrm,] #remove markers
    #print(sum(evidTab=="OL",na.rm=T))
    evidTab[evidTab=="OL"] = NA #ensuring no NA
    evidDat = sample_tableToList(evidTab)
    evidData = c(evidData,evidDat) #combine 
    
    repNames <- names(evidDat)
    if(length(repNames)>1) {
      for(rind in seq_along(repNames)) {
        sampleNameList[[length(sampleNameList)+1]] <- repNames[rind] 
      }
    }
    sampleNameList[[length(sampleNameList)+1]] <- repNames
  }
 return(list(popFreq=popFreq,evidData=evidData,refData=refData,sampleNameList=sampleNameList))
}
