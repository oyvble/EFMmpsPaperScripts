#script to compare degradation estimates obtained from different methods (compared against default)
#rm(list=ls())
calcfolds  = c("calcs","Supplementary")
methods0 = c("MAE1a","MAE1b","MAE2a","MAE2b")  #methods to consider
methods1 = paste0(methods0,"_ND")  #methods to consider
methodCompAgainst = "default" #comparison against this method
allMethods = c("Joint",methods0,methods1,methodCompAgainst) 

#Obtain relevant file names to import from
pat = "DCtables_"
filesImport = NULL
for(fold in calcfolds) {
  fns = list.files(fold,full=T,pattern=pat)
  filesImport = c(filesImport,fns)
}

#read estimates for each method:
par_list = list()
for(method in allMethods) {
#  method = allMethods[1]
  fns = basename(filesImport)
  idx = which(fns==paste0(pat,method,".RDS"))
  fn = filesImport[idx]
  
  #get DC res
  DClist <- readRDS(fn)
  sampleNames = names(DClist)
  par_list[[method]] = list()
  for(sampleName in sampleNames) {
  #  sampleName = sampleNames[1]
    DCobj = DClist[[sampleName]]
    noc = DCobj$NOC
    Mx = DCobj$theta[1:noc]
    th = DCobj$theta[-(1:noc)] #get estimates
    par_list[[method]][[sampleName]] = list(Mx=Mx,Th=th)
  }
} 

#Make comparison metric for each sample
sampleNames = unique(unlist(lapply(par_list,names)))
types = c("Degrad. slope","P.H.variability")
res_list = list() #
dumMat = matrix(NA,ncol=length(allMethods),nrow=length(sampleNames),dimnames = list(sampleNames,allMethods))
for(type in types) {
#  type=types[1]
  resMat = dumMat
  for(sampleName in sampleNames) {
    for(method in allMethods) {
      th = par_list[[method]][[sampleName]]$Th
      resMat[sampleName,method] =  th[grepl(type,names(th))]
    }
  }
  res_list[[type]] = resMat
}

#CREATING separate figures
m1 = "default" #compare against this
for(type in types) {
#  type = types[1]
  pdf(paste0(outfold,"//CompareParams_",type,"_against_",m1,".pdf"),height=6,width=6)
  for(m2 in setdiff(allMethods,m1)) {
  # m2 = allMethods[1] # "JointFast"
    v1 = res_list[[type]][,m1]
    v2 = res_list[[type]][,m2]
    
    #txt = paste0(m1," vs ",m2)
    txt = type
    xyrng = range(c(v1,v2))
    plot(v1,v2,main=txt,xlab = m1,ylab=m2,xlim=xyrng,ylim=xyrng);mtext(paste0("cor=",round(cor(v1,v2),2)));abline(a=0,b=1,lty=2)
    
    #Create colors etc
    isComb = grepl("/",names(v1))
    points(v1[isComb],v2[isComb],col=1,pch=19)
  }
  dev.off()
}

meths = allMethods[1]
for(i in 1:4) meths = c(meths,allMethods[c(i+1,i+5)])
#Put all in one plot
for (type in types) {
#  type = types[1]
  
  # ---- open ONE pdf page per type (one figure with many sub-figures) ----
  svg(file.path(outfold, paste0("CompareParams_", type, "_against_", m1, ".svg")),
      height = 10, width = 5)  # tweak size as you like
  
  # Layout: 5 rows x 2 cols; first plot spans both columns
  par(mfrow = c(5, 2),   # regular grid
      mar = c(4, 8, 1, 1),
      mgp =  c(2, 1, 0))
  
  # ---- draw panels ----
  for (i in seq_along(meths)) {
    m2 <- meths[i]
    v1 <- res_list[[type]][, m1]
    v2 <- res_list[[type]][, m2]
    xyrng <- range(c(v1, v2), na.rm = TRUE)
    
    plot(v1, v2, main = "", xlab = m1, ylab = m2, xlim = xyrng, ylim = xyrng)
    abline(a = 0, b = 1, lty = 2)
    isComb <- grepl("/", names(v1))
    points(v1[isComb], v2[isComb], col = 1, pch = 19)
    
    if(i==1 ) plot.new()   # empty panel
  }
  
  par(op)
  dev.off()
}

