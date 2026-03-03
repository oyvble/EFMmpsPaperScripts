
################
#help functions#
################

writeFile = function(x,fn,append=TRUE) write.table(t(x),fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=append)

getValres = function(x) { #store signif alleles
  if(sum(x$Signif)==0) return("")
  paste0(x[x$Signif,2],"-",x[x$Signif,3],"(",signif(x[x$Signif,5],4),")")
}  


importFreq = function(fn,format="RU",markignore=NULL) {
  sep="_"
  dat = readxl::read_xlsx(fn,sheet=1)
  cn = toupper(colnames(dat))
  locs = toupper(cn[-1]) #get loci names #Note to put space in "Penta X"
  av = unlist(dat[,1]) #get allele names

  popFreq = list()
  for(loc in locs) {
    #loc=locs[1]
    
    freq = unlist(dat[,cn==loc])
    induse = !is.na(freq)
    av2 = av[induse]
    freq2 = freq[induse]
    
    if(format != "LUS+") { 
      if(format=="RU") av2 = sapply( strsplit(av2,sep), function(x) x[1])
      if(format=="LUS") av2 = sapply( strsplit(av2,sep), function(x) paste0(x[1:2],collapse=sep))
    }
    agg = aggregate(freq2,by=list(av2),sum) #remember to aggregate RU also!
    av2 = agg[,1]
    freq2 = agg[,2]

    popFreq[[loc]] = freq2
    names(popFreq[[loc]]) = av2
  } #end for each loci
  if(!is.null(markignore)) popFreq = popFreq[setdiff(names(popFreq),markignore)] #keep those not in markignore
  return(popFreq)  
}


importEvid = function(fn,format="RU",markignore=NULL,threshT=NULL,skipEmpty=FALSE) {
  tab <- readxl::read_xlsx(fn, sheet = 1)
  cn = colnames(tab) #colnames 
  
  fixRU = function(x) { #helpfunction to avoid roundoff-errors on numbers
    suppressWarnings({ isnum = !is.na(as.numeric(x)) })
    x[isnum] = round(as.numeric(x[isnum]),1) 
    return(x)
  }
  
  sind = 2 #which(cn=="Renamed") #column with sample name 
  lind = 5 #which(cn=="Locus") #column with loci
  hind = 7 #which(cn=="Allele Coverage") #column with allele intensity
  samples <- unique(unlist(tab[,sind]))
  locs <- toupper(getKit("forenSeq","Marker")) #unique(unlist(tab[,lind]))
  if(!is.null(markignore)) locs = setdiff(locs,markignore)

  aindL = list(SEQ=11,RU=12,LUS=13,"LUS+"=14) #set index of different types
  aind = aindL[[format]] #get allele index  
 
#  dim(tab);length(Alist)
  ag <- aggregate(unlist(tab[,hind]),by=list(unlist(tab[,aind]) ,unlist(tab[,lind]),unlist(tab[,sind])),sum) #aggregate similar values (additivity)
  A  = ag[,1] #get alleles
  L = toupper(ag[,2]) #get loci
  S = ag[,3]#get sample names
  H = ag[,4] #get summed PHs
  
  datList = list() #structure samples into list (internal efm format)
  for(samp in samples) { #for each sample
#samp = samples[13]
    datList[[samp]] = list()
    for(loc in locs) { #for each locus
#loc = locs[23]  
      xind = which(S==samp & L==loc) #get index in X for given sample and locus
      if(length(xind)==0 && skipEmpty) next #DONT SKIP!
      keep <- which(!is.na(A[xind]) & A[xind]!="")

      PH <- as.numeric(as.character(H[xind][keep])) #get the peak heights
      if(!is.null(threshT)) keep = keep[PH>=threshT] #keep only alleles above thrshold (if given)

      datList[[samp]][[loc]]$hdata = PH[keep]
      datList[[samp]][[loc]]$adata = as.character(A[xind][keep])
    } #end for each loci
  } #end for each samples

  return(datList)
}

importRef = function(fn,format="RU",markignore=NULL) {
  tab <- readxl::read_xlsx(fn, sheet = 1)
  cn = colnames(tab) #colnames 
  
  fixRU = function(x) { #helpfunction to avoid roundoff-errors on numbers
    suppressWarnings({ isnum = !is.na(as.numeric(x)) })
    x[isnum] = round(as.numeric(x[isnum]),1) 
    return(x)
  }
  
  sind = which(cn=="Donor") #column with sample name 
  lind = which(cn=="Locus") #column with loci
  samples <- unique(unlist(tab[,sind]))
  locs <- toupper(getKit("forenSeq","Marker")) #unique(unlist(tab[,lind]))
  if(!is.null(markignore)) locs = setdiff(locs,markignore)
  
  renamedRefs = paste0("Ref",1:length(samples))
  aindL = list(SEQ=3,RU=4,LUS=5,"LUS+"=6) #set index of different types
  aind = aindL[[format]] #get allele index  
  Alist = unlist(tab[,aind])

  #  dim(tab);length(Alist)
  ag <- aggregate(unlist(Alist),by=list(Alist,unlist(tab[,lind]),unlist(tab[,sind])),length) #aggregate similar values (additivity)
  A  = ag[,1] #get alleles
  L = toupper(ag[,2]) #get loci
  S = ag[,3]#get sample names
  
  datList = list() #structure samples into list (internal efm format)
  for(samp in samples) { #for each sample
    datList[[samp]] = list()
    for(loc in locs) { #for each locus
#loc = locs[2]
#samp = samples[1]
      xind = which(S==samp & L==loc) #get index in X for given sample and locus
  #    if(length(xind)==0) next #DON*T SKIP!
      keep <- which(!is.na(A[xind]) & A[xind]!="")
      av = as.character(A[xind][keep])
      if(length(av)==1) av = rep(av,2)
      datList[[samp]][[loc]] = av
    } #end for each loci
  } #end for each samples
  names(datList) = renamedRefs #rename
  
  return(datList)
}

getData <- function(mixData2,refData2,popFreq) { #Helpfunction to get data to analyse
  locs <- names(popFreq)
  mixData <- lapply(mixData2,function(x) return(x[locs])) #return selected loci
  refData <- list()
#  for(loc in locs)  refData[[loc]] <- lapply(refData2,function(x) return(x[[loc]]$adata)) #return selected loci
  for(loc in locs)  refData[[loc]] <- lapply(refData2,function(x) return(x[[loc]])) #return selected loci
  
  #CALLING THE Qassignate function which converts non-observed alleles to "99" and includes missing allele frequencies
  Qret <- euroformix::Qassignate(samples=mixData, popFreq,  refData,incS=FALSE,incR=FALSE,minF=minF,normalize=FALSE )
  return(list(samples=mixData,refData=Qret$refData,popFreq=Qret$popFreq))
}

getVals = function(x) paste0( signif(x,3),collapse="/") #print helpfunction


#filtertype="static_withStutters";format="LUS+"
getAllData = function(filtertype="static_withStutters",format="LUS+", AT=30, skipEmpty=FALSE,datfold="data") {
  #  filtertypes = c("dynamic","static","static_withStutters")
#  formats = c("RU","LUS","LUS+")
  library(readxl)
  
  #Obtain references
  reffnFull = list.files(datfold,pattern="References.xlsx",full=T) #file with ref profiles
  popSel =  "AfAm" #select population (allele freqs)
  popfn = list.files(datfold,pattern=paste0("NIST 1036 ",popSel),full=T) #population freq file name
  
  evidfn = paste0("MixData30AT",filtertype,".xlsx") #file with evid profiles
  evidfnFull = list.files(datfold,pattern=evidfn,full=T)
  
  #READ ALL DATA AND RUN COMPARISONS:
  markignore = toupper(c("Amelogenin","Penta D")) #markers to not include
  evidData = importEvid(fn=evidfnFull,format,markignore,AT,skipEmpty) #mixData[[sample]][[locus]] = list(adata,hdata)
  refData = importRef(fn=reffnFull,format,markignore) #refData[[locus]][[sample]] = vector(allele1,allele2)

  #Include freqwuenc  
  popFreq = importFreq(fn=popfn,format,markignore) #read frequencies
  
  return(list(popFreq=popFreq,evidData=evidData,refData=refData))
}
