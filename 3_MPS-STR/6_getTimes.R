#script to analyze results:
#rm(list=ls())


#Obtain relevant file names to import from
method = "Joint"
pat = paste0("DCtables_",method)
fn = list.files(outfold,full=T,pattern=pat)

DCin = readRDS(file=fn)
samples = names(DCin) 
df = NULL
for(s in seq_along(samples)) {
#      s=1
  sample = samples[s]
  isComb = grepl("/",sample)
  time = DCin[[sample]]$timeMLE
  newrow = data.frame(Dataset=ds,Method=method,Sample=sample,isComb=isComb,Time=time)
  df = rbind(df,newrow)
}
df$NOC = substr(df$Sample,1,1)

#compareTimes(df,outfold) #use existing method
library(data.table)
setDT(df)
summary_dt <- df[, .(median_time = ceiling(median(Time, na.rm = TRUE))),  by = .(NOC,isComb)]
write.table(summary_dt,file=paste0(outfold,"//medianTime_",method,".csv"),sep=";",row.names = FALSE)

