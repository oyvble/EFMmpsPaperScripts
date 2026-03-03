
#type = "evid"

getAllData = function(AT=11,dynThresh=0.015,removeEmpty=TRUE, onlyREF=FALSE) {
#  removeEmpty=TRUE;AT=11;dynThresh=0.015
  getData = function(type) {
    library(readxl)
    evidfold = "Mixture (Evidence) Sample Reports"
    reffold = "Reference Sample Reports"
    sheets = c("Ancestry SNPs","Phenotype SNPs","Identity SNPs","Kinship SNPs")
    
    #evidence:
    fold = evidfold
    if(type=="ref") fold = reffold
    files = list.files(fold,full=TRUE)
    evidDat = list()
    for(f in seq_along(files)) {
      #  f=1
      file = files[f]
      evidName = strsplit(basename(file)," ")[[1]][1]
      evidDat[[evidName]] = list()
      for(sheet in sheets) {
        #      sheet = sheets[1]
        dat = read_excel(file,sheet,skip = 14)
        markers = unique(dat$Locus)
        for(marker in markers) {
          #        marker = markers[1]
          sdat = dat[dat$Locus==marker,]
          if(type=="ref") {
            geno = sdat$Genotype[1]
            if(geno=="INC") next #skip if not found
            adata = strsplit(geno,",")[[1]]
            evidDat[[evidName]][[toupper(marker)]] = adata
          } else {
            adata = sdat$`Allele Name`
            hdata = as.numeric(sdat$Reads)
            #          AT1 = ceiling(dynThresh*max(hdata))
            AT1 = floor(dynThresh*sum(hdata))
            indUse = hdata>=max(AT,AT1)
            if(sum(indUse)==0 && removeEmpty) next #skip marker if not considered
            if(sum(indUse)>2) indUse = order(hdata,decreasing = TRUE)[1:2] #keep 2 greatest
            evidDat[[evidName]][[toupper(marker)]] = list(adata=adata[indUse],hdata=hdata[indUse])
          }
          
        }
      }    
      #length(evidDat[[1]])
    }
    return(evidDat)
  }
  
  refData = getData("ref")
  if(onlyREF) return(list(refData=refData))
  
  
  #Read Data from reports
  library(CCMHr)
  popFreq = loadRDa("popFreq_1000G.rda")
  names(popFreq) = toupper(names(popFreq))
  locs_freqs <- names(popFreq)
  evidData = getData("evid")

  #check Number of alleles  
  #sapply(refData, function(x) max( sapply(x, function(y) length(unlist(y))) ) )
  #sapply(samplesAll, function(x) max( sapply(x, function(y) length(y$adata)) ) )
  
  #all(names(refData[[1]])==names(refData[[2]]))
  #REVERSE COMPLEMENT FREQUENCY TABLES IF NECESSARY
  reverseComp = setNames(c("A","T","C","G"),c("T","A","G","C"))
  for(loc in names(refData[[1]])) {
    #  loc = names(refData[[1]][1])
    refAlleles = refData[[1]][[loc]]
    freqAlleles = names(popFreq[[loc]])
    if(!all(refAlleles%in%freqAlleles)) { #Need reverse complement
      names(popFreq[[loc]]) = reverseComp[freqAlleles]
    }
    allAlleles  = unique(c(names(popFreq[[loc]]),
                           evidData[[1]][[loc]]$adata,refAlleles))
    if(length(allAlleles)!=2) stop("WRONG CONVERSION")
  }
  return(list(refData=refData,evidData=evidData,popFreq=popFreq))
}
