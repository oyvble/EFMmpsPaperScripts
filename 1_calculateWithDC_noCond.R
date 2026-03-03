#rm(list=ls())

#This is a main function for running the DC calculations (across all methods)
calculateWithDC_noCond = function(NOClist,sampleNameList,samplesAll,popFreq,modList,MAEoptions,methods,outfold = "calcs", MAEshrink=0,MAEshrink2=0) {
  dir.create(outfold,FALSE)
  
  #Obtain model variables (common for all samples)
  AT=modList$AT
  pC=modList$pC
  lambda=modList$lambda
  fst=modList$fst
  dynThresh = modList$dynThresh
  imputeDropout = modList$imputeDropout
  DEG=modList$DEG
  BWS=modList$BWS
  FWS=modList$FWS
  kit=modList$kit
  SNPmodule=modList$SNPmodule 
  nBins=modList$nBins
  
  for(mind in seq_along(methods)) {
    #mind=1
    method = methods[mind]
    print(method)
    DCfn <- paste0(outfold,"/DCtables_",method)
    MAEfile = paste0(outfold,"/MAEresults_",method) #saving MAE results to own file
    
    
    DClist <- list()
    markerEstimateList = list()
    for(sind in seq_along(sampleNameList)) { #for each samples in list
  #  sind = 28
      sampleNameVec = sampleNameList[[sind]]
      sampleNameReps <- paste0(sampleNameVec,collapse="/")
      print(sampleNameReps)
      #Obtain loci to be used:
      samplesUse = samplesAll[sampleNameVec]
      locs_sample <- unique(unlist(lapply(samplesUse,function(x) names(x))))
      if(!all(locs_sample%in%names(popFreq))) stop("Missing frequencies")
      NOC = NOClist[[sampleNameVec[1]]] #obtain NOC from defined list

      #Internal helpfunctions to do the MAE estimation
      helpEstimateMAE1 = function(useEmpirical=FALSE,inMAEshrink=0,forceNoDeg=FALSE) { #simple version based on sumPHs per marker
        kit2 = kit
        if(forceNoDeg) kit2 = NULL
        estimateMAE(samplesUse,kit=kit2,options=MAEoptions, imputeDropout=imputeDropout,useEmpirical=useEmpirical,seed=1, omega_max = 0.5, MAEshrink=inMAEshrink) #obtain estimates of MAE
      }
      
      helpEstimateMAE2 = function(useEmpirical=FALSE,inMAEshrink=0,inMAEshrink2=0,forceNoDeg=FALSE) { #simple version based on all alleles (and DC with MAEnaive)
        kit2 = kit
        if(forceNoDeg) kit2 = NULL
        estimateMAE2(NOC,samplesUse,popFreq,kit=kit2,DEG=DEG,BWS=BWS,FWS=FWS,AT=AT,pC=pC,lambda=lambda,fst=fst,dynThresh=dynThresh,
                     options = MAEoptions,imputeDropout=imputeDropout,useEmpirical=useEmpirical, seed=1, omega_max = 0.5, MAEshrink=inMAEshrink,MAEshrink2=inMAEshrink2 ) #obtain estimates of MAE
      }

      
      MAEmarker = NULL #no MAE estimate by default (also for binned)
      if(method=="binned") { #in case of using the binned approach (this is only for SNPs)
        mleHdNone = calcBinned(NOC,samplesUse,popFreq,AT=AT,pC=pC,lambda=lambda,fst=fst,dynThresh=dynThresh, doDC=TRUE,keepMLE=FALSE,SNPmodule=SNPmodule, nBins=nBins)
        DChdNone = mleHdNone$DC
        thetahat = mleHdNone$thetaMatrix #extract estimates
        Mx = colMeans(thetahat[,1:NOC,drop=FALSE])
        #  mleHdNone$thetaMatrix
        
      } else {         #Obtain estimates of MAE and calculate the time
        MAEtime = system.time({
          if(method=="MAE1a") MAEmarker = helpEstimateMAE1(TRUE) #ordinary MAE methods
          if(method=="MAE1b") MAEmarker = helpEstimateMAE1(FALSE) 
          if(method=="MAE2a") MAEmarker = helpEstimateMAE2(TRUE)
          if(method=="MAE2b") MAEmarker = helpEstimateMAE2(FALSE)
          if(method=="MAE1a_ND") MAEmarker = helpEstimateMAE1(TRUE,forceNoDeg=TRUE) #here we force no degradation info used when estimated MAE
          if(method=="MAE1b_ND") MAEmarker = helpEstimateMAE1(FALSE,forceNoDeg=TRUE) 
          if(method=="MAE2a_ND") MAEmarker = helpEstimateMAE2(TRUE,forceNoDeg=TRUE)
          if(method=="MAE2b_ND") MAEmarker = helpEstimateMAE2(FALSE,forceNoDeg=TRUE)
          if(method=="MAE1a+") MAEmarker = helpEstimateMAE1(TRUE,MAEshrink) #MAE methods with shrinkage (only empirical methods)
          if(method=="MAE2a+") MAEmarker = helpEstimateMAE2(TRUE,MAEshrink,MAEshrink2)
        })[3]
        if(!is.null(MAEmarker)) {
          markerEstimateList[[sampleNameReps]] = list(Est=MAEmarker,Time=MAEtime)
          saveRDS(markerEstimateList,file=paste0(MAEfile,".RDS"))
        }
        
        time = system.time({
          if(method=="Joint") { #This is joint inference of theta and MAE
            mleHdNone = calcJointMAP_full(NOC,samplesUse,popFreq,kit=kit, DEG=DEG, BWS=BWS, FWS=FWS, AT=AT,pC=pC,lambda=lambda,fst=fst,dynThresh=dynThresh, steptol=1e-4, hyppar_tau=c(2,1),hyppar_sigma=c(0,0.5,1))#$model$MAEmarker
          } else {
            #Fit model given estimated MAE and Deconvolute (without any conditioning)
            mleHdNone = calcMLE(NOC,samplesUse,popFreq,NULL,NULL, NULL,kit, DEG, BWS, FWS, AT,pC,lambda,fst,dynThresh=dynThresh,MAEmarker=MAEmarker,SNPmodule=SNPmodule)#,hyppar_sigma = c(0,0.5,1))
          }
        })[3]
        thetahat = mleHdNone$fit$thetahat2 #extract estimates
        if(any(is.na(thetahat))) next #skip deconvolution if not valid
        DChdNone = deconvolve(mleHdNone,maxlist = 1000,alpha=1) #results not knowing POI reference
        Mx=thetahat[1:NOC]
      }
      
      #STORING DC RESULT
      DClist[[sampleNameReps]] = list(MargGeno=DChdNone$table3, MargAllele=DChdNone$table4,
                                          NOC=NOC,theta=thetahat,Mx=Mx, condRef = NULL, timeMLE=time)
      saveRDS(DClist,file=paste0(DCfn,".RDS"))
    } #end for each evidence
  } #end for each method
}


######################################################################





