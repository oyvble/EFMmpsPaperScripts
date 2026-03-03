#THIS SCRIPT WILL REPRODUCE RESULTS FOUND IN PAPER/SUPPLEMENTARY
#rm(list=ls())
#library(EFMmps)
maindir = "~\\EFMmpsPaperScripts" 
setwd(maindir)
source("2_compareWithTrueRefs.R") #makeGrouping function used to compare MAE classifications
source("4_createSummaryResults.R") #helpfunction for comparing results
source("5_createOtherResults.R") #helpfunction for creating other results
outfold = "paperResults"
calcfold = "calcs"
resfold = "Results"
dir.create(outfold,showWarnings = F)

#Obtain folders with the different datasets
dsets = list.dirs(recursive = FALSE,full.names = FALSE)#pattern="_")
dsets = dsets[grepl("_",dsets)] #
MAEmethods = c("MAE1a", "MAE1b", "MAE2a","MAE2b")
methods = c(MAEmethods,"default","binned")
datNames = sapply(strsplit(dsets,"_"),function(x) x[2])

df = NULL
for(d in seq_along(dsets)) {
  #  d=4
  ds = dsets[d]
  datName = datNames[d]
  dfin = readRDS(file=paste0(ds,"/",resfold,"/df_SummaryResults.RDS"))
  newdf = cbind(Dataset=datName,dfin)
  df = rbind(df,newdf)
}
#View(df)
df$isCombined[df$Dataset=="MH"] = FALSE
dfCE = df[df$Dataset=="CE" & df$isCombined=="All",]
dfCE$isCombined = TRUE
dfNotCE = df[df$Dataset!="CE",]
dfNotCE = dfNotCE[dfNotCE$isCombined!="All",]
df = rbind(dfNotCE,dfCE)
#df0 = df #copy
#df = df0
#View(df)
df$Dataset = factor(df$Dataset,levels=datNames)
df$Method = factor(df$Method,methods)
df = df[!is.na(df$Method),] #ensure that no other methods are included
df$isCombined = factor(df$isCombined,c(TRUE,FALSE))
df$combined <- factor(df$isCombined, labels = c("Combined","Non-combined"))
#Include true ref and compare with DC results

metrics = c("Brier","Accuracy","Coverage","Acc_high")
for(met in metrics) {
  makeSummaryFig(df,met,outfold)
  
  sdf = df[df$Group != "Hard",]
  makeSummaryFig(sdf,met,outfold,"notHard")
}


#OBTAIN MAE COMPARISON RESULTS
df <- NULL
for(d in seq_along(dsets)) {
  #  d=4
  ds = dsets[d]
  datName = datNames[d]
  for(m in seq_along(MAEmethods)) {
# m=1
    method = MAEmethods[m]
    pat = paste0("MAEresults_",method)
    fn = list.files(paste0(ds,"/",calcfold),pattern=pat,full=T)
    dfin = readRDS(fn) #read results
    
    sampleNames = names(dfin) #get calculated results
    for(sampleName in sampleNames) {
#      sampleName=sampleNames[1]
      MAEhat = dfin[[sampleName]]$Est
      MAEtime = dfin[[sampleName]]$Time
      locs = names(MAEhat) 
      est = as.numeric(MAEhat)
      isComb = grepl("/",sampleName)
      if(datName=="CE") isComb=FALSE #force 
      newdf = data.frame(Dataset=datName,Method=method,Sample=sampleName,isComb=isComb,Locs=paste0(locs,collapse="/"),MAE=paste0(est,collapse="/"),Time=MAEtime)
      df = rbind(df,newdf)
    }
  }
}
compareMAEresults(df, outfold) #create plots that compare MAE between methods
compareTimes(df, outfold) #create plots that compare Time between methods


#View(df)

#########################################################
##ANALYSE HOW WELL MAE CLASSIFIES THE DIFFICULTY GROUPS##
#########################################################

#Obtain folders with the different datasets
#dsets = list.dirs(recursive = FALSE,full.names = FALSE)#pattern="_")
#dsets = dsets[grepl("_",dsets)] #
dsets = dsets[!grepl("CE",dsets)] #dont include CE
datNames = sapply(strsplit(dsets,"_"),function(x) x[2])
methods = c("MAE1a", "MAE1b", "MAE2a","MAE2b")
levels = c("Easy","Medium","Hard")

df = NULL #obtain comparison table (true based vs estimated based)
for(d in seq_along(dsets)) {
  #  d=4
  ds = dsets[d]
  datName = datNames[d]
  
  for(method in methods) {
    pattern = paste0("DCtables_",method)
    dcFile = paste0(ds,"/",calcfold,"/",pattern,".RDS")
    resList = readRDS(file=paste0(ds,"/",calcfold,"/estimatedMxWithTrue.RDS"))   #Dont restrict to combined rule
    groupList = readRDS(file=paste0(ds,"/",calcfold,"/ruleList.RDS"))   #Dont restrict to combined rule
    
    DClist = readRDS(file=dcFile)  
    sampleNames = names(DClist) 
    for(s in seq_along(sampleNames)) {
      #      s=2
      sampleName = sampleNames[s]
      MxHat = DClist[[sampleName]]$Mx #obtaine estimated based on MAE
      MxTrue = resList[[sampleName]][,"Mx"]
      NOC = length(MxTrue)
      if(is.null(MxTrue)) next
      MxTrueSortedInd = order(resList[[sampleName]][,"Mx"])
      resList[[sampleName]][MxTrueSortedInd,"Mx"] = sort(MxHat) #insert same sorting
      maxDiff = max(abs(sort(MxTrue) - sort(MxHat)))
      
      #Check if same group was obtained
      updatedGroup = makeGrouping(resList[sampleName])[[1]]
      ordDonors1 = match(names(updatedGroup),rownames(resList[[sampleName]]))
      ordDonors2 = match(names(groupList[[sampleName]]),rownames(resList[[sampleName]]))
      groupMAE = factor(updatedGroup[ordDonors1],levels=levels,ordered=TRUE)
      groupTrue = factor(groupList[[sampleName]][ordDonors2],levels=levels,ordered=TRUE)
      groupDiff = as.integer(groupMAE) - as.integer(groupTrue)
      groupDiffComb = paste0(groupDiff,collapse="/")
      groupDiffTot = sum(groupDiff!=0)
      groupDiffnSerious = sum(groupDiff<0) #Pred as easier than with true is serious
      
      #groupTrueTxt = paste0(groupTrue,collapse="/")
      #groupMAETxt = paste0(groupMAE,collapse="/")
      groupTrueTxt = paste0(substr(groupTrue,1,1),collapse="/")
      groupMAETxt = paste0(substr(groupMAE,1,1),collapse="/")
      MxTrueTxt = paste0(round(MxTrue,2),collapse="/")
      MxHatTxt = paste0(round(resList[[sampleName]][,"Mx"],2),collapse="/")
      newrow = data.frame(Dataset=datName,Method=method,Sample=sampleName,isComb=grepl("/",sampleName),mxTrue=MxTrueTxt,mxHat=MxHatTxt,mxMaxDiff=maxDiff,groupMAE=groupMAETxt,groupTrue=groupTrueTxt,groupDiff=groupDiffComb,groupDiffTot=groupDiffTot,groupDiffnSerious=groupDiffnSerious,NOC=NOC)
      df = rbind(df,newrow)
    }
  }
}
#View(df)
library(dplyr)

#0) total number of predictions
sum(df$NOC[df$Method=="MAE1a"]) #292


#1) SUMMARIZE TOTAL NUMBER OF DIFFERENCS PER METHOD
df_summary_method <- df %>%
  group_by(Method) %>%
  summarise(totalSerious = sum(groupDiffnSerious, na.rm = TRUE),totalDiffs = sum(groupDiffTot, na.rm = TRUE), N=n(),.groups = "drop")

#Store table
write.table(df_summary_method,file=paste0(outfold,"/predictDifficultyGroup_compareMAE.csv"),sep=";",row.names = F)

#2) SHOW MORE DETAILS WHAT DATASETS THESE DIFFERENCES ARE FOR
df_summaryList = list()
for(method in methods) {
  df_summary <- df %>%filter(Method==method)%>%
    group_by(Dataset, isComb) %>%
    summarise(totalSerious = sum(groupDiffnSerious),totalDiffs = sum(groupDiffTot, na.rm = TRUE), N=n(),.groups = "drop")
  df_summaryList[[method ]] = df_summary[df_summary$totalDiffs>0,]
  
}
View(df_summaryList$MAE2b)

#3) Print tables with more details:
showCols = c("Dataset","Sample","mxTrue","mxHat","groupTrue","groupMAE")
names(df_detail)
for(method in methods) {
  
  df_detail <- df %>%filter(Method==method, groupDiffTot>0)%>%group_by(Dataset, Sample)
  out = df_detail[,showCols]
  
  #Store table
  write.table(out,file=paste0(outfold,"/predictDifficultyGroup_details_",method,".csv"),sep=";",row.names = F)
}


