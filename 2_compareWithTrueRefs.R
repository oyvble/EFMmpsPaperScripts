
#This is the rule based
makeGrouping = function(MxResList,restrictAsCombined=FALSE) {
  #Create rules first: Obtain Easy, Middle, Hard
  dropThresh = c(0.005,0.1) #accept some dropouts 0.5% (necessary for SNPs)
  levels = c("Easy","Medium","Hard")
  rule_list = list() #include details about samples
  for(sampleName in names(MxResList)) {
    #  sampleName=names(MxResList)[2]
    MxResList_sample = MxResList[[sampleName]]
    MxResList_sample = MxResList_sample[order(MxResList_sample[,1],decreasing = TRUE),,drop=F] #sort by Mx
    refNames = rownames(MxResList_sample)
    NOC = length(refNames)
    Mxhat_true = MxResList_sample[,1] #Mx based on true refs
    dropProb_true = MxResList_sample[,2] #dropout rate of true ref
    category = rep(3,NOC) #All hard by default
    
    #1st component:
    mxDiff1 = Mxhat_true[1] - Mxhat_true[2] #take difference between largest and second largest
    if(mxDiff1>=0.2) {
      category[1] = 1 #set as Easy
    } else if(mxDiff1>=0.1) {
      category[1] = 2 #set as Medium
    }
    
    #Additional requirement for 1st component when NOC>2
    if(NOC>2 && Mxhat_true[1]<0.6) { #if 1st below 60%
      category[1] = ifelse(Mxhat_true[1]>=0.5,2,3) #set to medium or hard if over/under 50%
    } 
    
    #2nd component:
    if(NOC==2) {
      category[2] = 1 #set same category as 1st (depended on Mx)
      
    } else { #otherwise we have 3p mixture (check distance with Mx3)
      mxDiff2 = Mxhat_true[2] - Mxhat_true[3] #take difference between second and third largest component
      if(mxDiff2>=0.2) {
        category[2] = 1 #set as Easy
      } else if(mxDiff2>=0.1) {
        category[2] = 2 #set as Medium
      }
      
      #3rd component:
      category[3] = category[2] #set same category as 2nd (depended on Mx)
    }
    

    #Dropout requirement
    indDrop1 = dropProb_true>=dropThresh[1]
    indDrop2 = dropProb_true>=dropThresh[2]
    category[indDrop1] = pmax(category[indDrop1], 2)  #set as medium
    category[indDrop2] = pmax(category[indDrop2], 3)  #set as Hard
    category[Mxhat_true<0.1] = 3 #Mx below 10% always as Hard
    
    #Decreasing group requirement
    #Loop and set maxium level as all previous components
    for(c in 2:NOC) {
      category[c] = max(category[c], category[1:(c-1)])
    }
    categoryName = levels[category]
    rule_list[[sampleName]] = setNames(categoryName,refNames)
  }
  
  if(restrictAsCombined) {
    sampleNames = names(rule_list)
    isCombined = grepl("/",sampleNames)
    if(any(isCombined)) {
      samplesCombined = sampleNames[isCombined]
      for(s in seq_along(samplesCombined)) {
        sampleList = strsplit(samplesCombined[s],"/")[[1]]
        for(r in seq_along(sampleList)) {
          rule_list[[sampleList[[r]]]] =  rule_list[[samplesCombined[s]]] #set as the same
        }        
      }
    }
  }
  
  return(rule_list)
}

#Helpfunction that estimate Mx based on normalized signal (simple OLS model)
#estimateMxWithTrueRefs(sampleNameList,evidList,refList,trueRefList,calcFold,AT=50,dynThresh=0)
#evidList=samplesAll;
estimateMxWithTrueRefs = function(sampleNameList,evidList,refList,trueRefList,calcFold,AT,dynThresh=0) {
  imputeDropout = AT/2 #impute read if missing for true contr.
   
  resList <- list()
  for(sind in seq_along(sampleNameList)) { #for each samples in list
#  sind = 6
    sampleNameVec = sampleNameList[[sind]]
    sampleNameReps <- paste0(sampleNameVec,collapse="/")
    #Obtain loci to be used:
    samplesUse = evidList[sampleNameVec]
    locs_sample <- unique(unlist(lapply(samplesUse,function(x) names(x))))
    
    #Calculate while condition on all true
    trueRefs = trueRefList[[sampleNameVec[1]]]
    refDatTrue = refList[trueRefs]
    NOC = length(refDatTrue) #length(strsplit(sampleNameVec[1],"_")[[1]]) #Obtain number of contributor based on name

    locs_use = locs_sample
    for(r in seq_along(refDatTrue)) {
      locs_use = intersect(locs_use,names(refDatTrue[[r]]))
    }
#    locs_refs = unique(unlist(lapply(refDatTrue, function(x) names(x))))
 #   locs_use = intersect(locs_sample,locs_refs)
    length(locs_use)
    
    #Part 1: Obtain LCB and Contribution overlap probabilities
    #Obtain elements to calculate linear model fit
    propShared = matrix(NA,nrow=length(locs_use),ncol=NOC) #calculated for each marker and contributor
    nDropout = matrix(0,nrow=length(locs_use),ncol=NOC) #calculated for each marker and contributor
    Yresponse = NULL
    Xmat = NULL
    for(l in 1:length(locs_use)) {
#      l=1
      loc=locs_use[l]
      refLocs = lapply(refDatTrue,function(x) unlist(x[[loc]]))
      if(any(sapply(refLocs,is.null))) next #skip marker if not all refs given
      evidLocs = lapply(samplesUse,function(x) x[[loc]])
      nR = length(evidLocs) #number of replicates
      if(is.null(dynThresh)) {
        evidAlleles = lapply(evidLocs,function(x) x$adata) #already imputed
        evidHeights = lapply(evidLocs,function(x) x$hdata)
      } else {
        if(dynThresh>0) {
          stop("NOT IMPLEMENTED")
        } else {
          evidAlleles = lapply(evidLocs,function(x) x$adata[x$hdata>=AT])
          evidHeights = lapply(evidLocs,function(x) x$hdata[x$hdata>=AT])
        }
      }
      uniqueAllelesEvids = unique(unlist(evidAlleles))
      uniqueAllelesRefs = unique(unlist(refLocs))
      uniqueAlleles = sort(uniqueAllelesRefs)#sort(intersect(uniqueAllelesEvids,uniqueAllelesRefs))
      
      contrMat = NULL    #establish true contribution   
      for(c in 1:length(refLocs)) {
        #        c=2
        contrMat = cbind(contrMat,table( factor(refLocs[[c]],levels=uniqueAlleles)))
        
        #Calculate proportion of shared alleles
        allelesC = refLocs[[c]]
        isNotDropout = allelesC%in%uniqueAllelesEvids #indicate if allele of ref is observed
        nDropout[l,c] = sum(!isNotDropout)
        allelesC = allelesC[isNotDropout] #avoid calculate on missing alleles
        allelesOthers = unlist(refLocs[-c])
        sharedAlleles = sum(allelesC%in%allelesOthers)
        totalAlleles = length(allelesC)
        propShared[l,c] = sharedAlleles/totalAlleles
      }
      #colnames(contrMat) = paste0("C",1:ncol(contrMat))
      contrMat = contrMat/2 #adjust
      
      for(r in 1:nR) {
        #        r=1
        yv = setNames(rep(0,length(uniqueAlleles)),uniqueAlleles)
        isNotMissing = evidAlleles[[r]]%in%uniqueAlleles
        yv[!uniqueAlleles%in%evidAlleles[[r]]] = imputeDropout
        yv[evidAlleles[[r]][isNotMissing]] = evidHeights[[r]][isNotMissing]
        #     yv = log(yv)
        yv = yv/sum(yv) #normalize
        Yresponse = c(Yresponse,yv)
        Xmat = rbind(Xmat,contrMat)
      }
    } #end for each locus
    #  dim(Xmat)
    #length(Yresponse)
    propSharedMeans = colMeans(propShared,na.rm = T) #get proportion of shared alleles (per contributor)
    propDropMeans = colMeans(nDropout,na.rm = T) #get proportion of shared alleles (per contributor)
    propDropMeans = propDropMeans/2 #number of outcome is 2 per marker (must be adjusted for)
    
    #estimate Mx and sigma
    fit = summary(lm(Yresponse~Xmat-1)) #fit without intercept
    MxHat = fit$coefficients[,1]
    MxHat = MxHat/sum(MxHat) #normalize
    sigmaHat = fit$sigma
    
    
    score = cbind(Mx=MxHat,propDropMeans=propDropMeans,propSharedMeans=propSharedMeans)
    rownames(score) = trueRefs
    resList[[sampleNameReps]] = score

  } #end for each sample
  saveRDS(resList,file=paste0(calcFold,"/estimatedMxWithTrue.RDS"))
  saveRDS(makeGrouping(MxResList=resList) ,file=paste0(calcFold,"/ruleList.RDS"))   #Dont restrict to combined rule
  saveRDS(makeGrouping(MxResList=resList,restrictAsCombined=TRUE),file=paste0(calcFold,"/ruleList_asCombined.RDS")) #restrict
  
}

#Make a table of Mx and PHvar results from DC results  
makeMxResults = function(calcFold = "calcs",outFold="Results") {
  library(ggplot2)
  #Include DC results
  fnsFull = list.files(calcFold,full=T,pattern="DCtables")
  methods = sapply(strsplit(basename(fnsFull),"_"), function(x) gsub(".RDS","",x[2]))
  
  helper = function(x) paste0(sort(x,decreasing = TRUE),collapse="/")
  
  dfPar = NULL
  for(m in seq_along(methods)) {
    # m=1
    method=methods[m]
    #   print(method)
    fn = fnsFull[m]
    
    #get DC res
    DClist <- readRDS(fn)
    sampleNames = names(DClist)
    for(sampleName in sampleNames) {
      #    sampleName=sampleNames[1]
      
      #    sampleName = "3P_0.067ng_10-2-1"
      DCobj = DClist[[sampleName]]#$NONE
      if(method=="binned") {
        thetaMat=DCobj$theta
        MxRange = grepl("Mix-prop",colnames(thetaMat))
        Mx = colMeans(thetaMat[,MxRange])
        Mx0 = Mx/sum(Mx)
        Mx = helper(round(Mx0,2))
        PHexp = round(mean(thetaMat[,max(which(MxRange)) + 1]),3)
        PHvar = round(mean(thetaMat[,max(which(MxRange)) + 2]),3)
      } else {
        Mx0 = DCobj$Mx
        Mx = helper(round(Mx0,2))
        PHexp = round(DCobj$theta["P.H.expectation"],3)
        PHvar = round(DCobj$theta["P.H.variability"],3)
      }
      PHcontr = helper(round(PHexp*Mx0,1))
      newrows = data.frame(Method=method,Sample=sampleName,Mx=Mx,PHvar=PHvar,PHcontr=PHcontr)
      dfPar = rbind(dfPar,newrows)
    } #end each sample
  } #end for each method
  rownames(dfPar) = NULL
  saveRDS(dfPar,file=paste0(calcFold,"/ParamResults.RDS"))
  #return(df)
  
  #Study parameters
  #dfPar = df #makeMxResults()
  rename <- renameHead <- sort(unique(dfPar$Sample))
  isCombined = grepl("/",rename)
  rename[isCombined] = sapply( strsplit(rename[isCombined],"/"),function(x) substr(x[1],1,nchar(x[1])-1))
  names(rename) = renameHead
  
  renameNoComb = rename[!isCombined]
  renameNoComb = renameNoComb[order(nchar(renameNoComb))]
  renameComb = rename[isCombined]
  renameComb = renameComb[order(nchar(renameComb))]
  
  renameSorted = c(renameNoComb,renameComb)
  dfPar$SampleNew = factor(rename[dfPar$Sample],levels=renameSorted)
  methodLevels = c("default","MAE1a","MAE1b","MAE2a","MAE2b" ,"binned")
  #if(!all(unique(dfPar$Method)%in%methodLevels)) stop("A method was not found")
  dfPar = dfPar[dfPar$Method%in%methodLevels,] #make subset of data
  methodLevelsUse = methodLevels[methodLevels%in%dfPar$Method]
  dfPar$Method = factor(dfPar$Method,levels=methodLevelsUse)
  
  #pdf(paste0(outFold,"/Figure_PHvarResults.pdf"),width=25,height=9)
  gg = ggplot(dfPar, aes(x = Method, y = PHvar, fill = Method)) +  geom_boxplot() +
    labs(x = "Method", y = "PHvar", title = "Boxplot of PHvar by Method") + theme(legend.position = "none")
  ggsave(filename = paste0(outFold,"/Figure_PHvarBarplot.svg"),plot=gg,width=6,height=6)
  
  gg = ggplot(dfPar, aes(x = SampleNew, y = PHvar, color = Method, group = Method)) +
    geom_line() + geom_point() + labs(x = "Sample", y = "PHvar", title = "Trend of PHvar per Sample by Method") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
  ggsave(filename = paste0(outFold,"/Figure_PHvarTrendplot.svg"),plot=gg,width=12,height=6)
  
  
  Mx = dfPar$Mx
  dfPar$Mx1 = sapply(strsplit(Mx,"/"),function(x) paste0(round(as.numeric(x),1),collapse="/"))
  
  gg2 = gg + geom_text(aes(label = dfPar$Mx1),  vjust = -0.5, size = 2,color=1)        # adjust text size
  ggsave(filename = paste0(outFold,"/Figure_PHvarTrendplotWithMx.svg"),plot=gg2,width=12,height=6)
  # return(dfPar)
} #done function


compareDCWithTrueRefs = function(refList,calcFold = "calcs") {
  library(data.table)
  MxResList = readRDS(file=paste0(calcFold,"/estimatedMxWithTrue.RDS"))
  
  #Include DC results
  fnsFull = list.files(calcFold,full=T,pattern="DCtables")
  methods = gsub("DCtables_","",basename(fnsFull)) #sapply(strsplit(basename(fnsFull),"_"), function(x) gsub(".RDS","",x[2]))
  methods = gsub(".RDS","",methods)
  helper = function(x) paste0(sort(x),collapse="/")
  
  summary_df = list() #obtain summary table per sample (estimated Mx, group, Top Accuracy)
  res_list = list() #this is DC based results (match and probs)
  for(m in seq_along(methods)) {
# m=5
    method=methods[m]
    print(method)
    fn = fnsFull[m]
  
    #get DC res
    DClist <- readRDS(fn)
    sampleNames = names(DClist)
    res_list[[method]] = list() #for storing results
    for(sampleName in sampleNames) {
  #    sampleName=sampleNames[1]
      res_list[[method]][[sampleName]] = list()
      
  #    sampleName = "3P_0.067ng_10-2-1"
      DCobj = DClist[[sampleName]]#$NONE
      
      #Obtain Mx estimates:
      MxResList_sample = MxResList[[sampleName]]
      if(is.null(MxResList_sample)) next #skip marker if not to evaluate
      MxResList_sample = MxResList_sample[order(MxResList_sample[,1]),,drop=F] #sort by Mx
      Mxhat_true = MxResList_sample[,1] #Mx based on true refs
      Mxhat_DC = DCobj$Mx #Mx based on model
      
    #  View(DCobj$MargGeno)
      C_IDs = unique(DCobj$MargGeno[,1])
      C_IDs = C_IDs[C_IDs!=""]
      NOC = length(C_IDs)
      if(method=="binned") Mxhat_DC = colMeans(DCobj$theta[,1:NOC,drop=FALSE])
      order_DC = order(Mxhat_DC) #also get order of DC
      
      #Create a matchmatrix between predictions and refs    
      tabDC = DCobj$MargGeno
  #    head(tabDC)

      refNames = names(Mxhat_true)
      probList = list() #matrix(list(),ncol=1,nrow=length(refNames),dimnames = list(refNames,NULL))
      matchVec = rep(0,NOC)
      locusCounterVec = rep(0,NOC)
      for(r in seq_along(refNames)) { #traverse each
#        r=1
        refName = refNames[r]
        refDat = refList[[refName]]
        locs = names(refDat)
        cUse = order_DC[r] #obtain DC contr ID to use (based on model based Mx, ordered)

        C_ID=C_IDs[cUse]
        probMat = NULL
        for(loc in locs) {
      #  loc=locs[2]
          trueGeno = unlist(refDat[[loc]])#$adata
          if(length(trueGeno)==0) next #skip if no data
          if(length(trueGeno)==1) trueGeno = rep(trueGeno,2) #give as twice if homozygous
          
          #obtain predicted genotypes          
          #tabDC[tabDC[,2]==loc,]
          rows = tabDC[,1]==C_ID & tabDC[,2]==loc
          if(!any(rows)) next #skip if no DC result found
          subtabDC = tabDC[rows,,drop=FALSE]
          allAllelesList = strsplit(subtabDC[,3],"/")
          allAlleles = unique(unlist(allAllelesList))
          if(all(allAlleles=="99")) next #skip if all DC alleles are dropout (no data)
          
          #Continue:
          trueGeno[!trueGeno%in%allAlleles] = "99" #assign as dropout
          #CHECK IF ALLELES ARE WITHIN DATA
          trueGeno = helper(trueGeno)#get true genotype
          predGenos = sapply( allAllelesList, helper)          

          #UPDATE COUNTING            
          matchVec[r] = matchVec[r] + as.integer(predGenos[1]==trueGeno) #whether top is a match
          locusCounterVec[r] = locusCounterVec[r] + 1
          
              #indicicate which is matching
          indTrueGeno = predGenos==trueGeno
          probNewRows = data.frame(Ref=refName,Locus=loc,Prob=as.numeric(subtabDC[,4]), Match=0)
          probNewRows$Match[indTrueGeno] = 1  
          probMat = rbind(probMat,probNewRows)
        } #end each locus
        probList[[refName]] = probMat
      } #end each ref
      
      res_list[[method]][[sampleName]] = probList
      
      Accuracy = round(matchVec/locusCounterVec,3)
      MxVec = round(MxResList[[sampleName]][refNames,"Mx"],3)
      DropVec = round(MxResList[[sampleName]][refNames,"propDropMeans"],3)
      newrows = data.frame(Sample=sampleName, Ref=refNames,Mx=MxVec,Drop= DropVec,Method=method,Acc=Accuracy,nMatches=matchVec,nMarkers=locusCounterVec)
      summary_df = rbind(summary_df,newrows)
    } #end each sample
  } #end for each method
  rownames(summary_df) = NULL
  saveRDS(summary_df,file=paste0(calcFold,"/ComparedWithTrue_summaryTable.RDS"))
  
  # 1) Bind by ref within each (method, sample)
  L1 <- lapply(res_list, function(level_sample) {
    lapply(level_sample, function(level_ref) {
      rbindlist(level_ref, idcol = "Ref", use.names = TRUE, fill = TRUE)
    })
  })
  
  # 2) Bind by sample within each method
  L2 <- lapply(L1, function(level_sample_bound) {
    rbindlist(level_sample_bound, idcol = "Sample", use.names = TRUE, fill = TRUE)
  })
  
  # 3) Bind by method
  df <- rbindlist(L2, idcol = "Method", use.names = TRUE, fill = TRUE)
#View(df)
  saveRDS(df,file=paste0(calcFold,"/ComparedWithTrue_allProbs.RDS"))
  #return(df)
} #done function

