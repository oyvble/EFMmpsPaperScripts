
compareMAEresults = function(df, outfold) {
  
  ds = unique(df$Dataset)
  metNames = unique(df$Method)
  #ms = unique(df$Method)
  
  corrList = list()
  slopeList = list()
  dfList = list()
  for(d in ds) {
  #  d=ds[4]
    print(d)
    sdf = df[df$Dataset==d,]
    samples = unique(sdf$Sample)
    corrList[[d]] = list()
    slopeMat = NULL
    dfPlot = NULL
    for(sample in samples) {
    #  sample=samples[10]
      ssdf = sdf[sdf$Sample==sample,]
      vals = NULL
      for(method in metNames) {
  #      method = metNames[3]
        sub = ssdf[ssdf$Method==method,]
        locs = strsplit(sub$Locs,"/")[[1]]
        mae = as.numeric(strsplit(sub$MAE,"/")[[1]])
        vals = cbind(vals,setNames(mae,locs))
        
        dfPlot = rbind(dfPlot,data.frame(Dataset=d,Sample=sample,Method=method,Marker=locs,MAE=mae))
      }
      colnames(vals) = metNames
      corrList[[d]][[sample]] =  cor(vals)
      f12 = lm(vals[,2]~vals[,1])
      f34 = lm(vals[,4]~vals[,3])
      slopeMat = rbind(slopeMat,c(f12$coefficients[2],f34$coefficients[2]))
      #    plot(vals[,1],vals[,2]);abline(a=0,b=1)
  # plot(vals[,3],vals[,4]);abline(a=0,b=1)
    }
    dfList[[d]] = dfPlot
    rownames(slopeMat) = samples
    colnames(slopeMat) = c("slope12","slope34")
    slopeList[[d]] = slopeMat
    
  }
  library(data.table)
  mae_long <- rbindlist(dfList)#, idcol = "Dataset")
  #head(dfPlot)
  
  corrRange = NULL
  slopeRange = NULL
  for(d in ds) {
  #  d=ds[1]
    corrElem = corrList[[d]]
    corr12 = sapply(corrElem,function(mat) mat[1,2])
    corr34 = sapply(corrElem,function(mat) mat[3,4])
    
    makeTxt = function(x) {
      x = round(x,3)
      paste0("[",x[1]," - ",x[2],"]")
    } 
    row = c(makeTxt(x=range(corr12,na.rm = TRUE)),makeTxt(x=range(corr34,na.rm = TRUE)))
    corrRange = rbind(corrRange,row)
  
    slopeElem = slopeList[[d]]
    row = c(makeTxt(x=range(slopeElem[,1],na.rm = TRUE)),makeTxt(x=range(slopeElem[,2],na.rm = TRUE)))
    row2 = round(colMeans(slopeElem),2)
    slopeRange = rbind(slopeRange,c(row,row2))
  }
  rownames(corrRange) <- rownames(slopeRange) <- ds
  colnames(corrRange) <- c("MAE1a vs MAE1b","MAE2a vs MAE2b")
  colnames(slopeRange) <- c("MAE1a vs MAE1b","MAE2a vs MAE2b","mean1","mean2")
  
  write.table(corrRange,file=paste0(outfold,"//MAEcompare_corRange.csv"),sep=";")
  write.table(slopeRange,file=paste0(outfold,"//MAEcompare_slopeRange.csv"),sep=";")
  
  #Show scatterplot relation for each datasets:
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # assume your object is called `mae_long`
  # columns: Dataset, Sample, Method, Marker, MAE
  #mae_long_small <- mae_long
  # 1) Decide how many to sample per Dataset/Sample
  set.seed(123)
  mae_long_small <- mae_long %>%
    group_by(Dataset, Sample) %>%
    group_modify(function(df, key) {
      n_row <- nrow(df)
      # sample up to 1000 rows per group
      df[sample.int(n_row, size = min(5000L, n_row)), ]
    }) %>%
    ungroup()
  
  # step 1: wide format per Dataset–Sample–Marker
  mae_wide <- mae_long_small %>%
    select(Dataset, Sample, Marker, Method, MAE) %>%
    pivot_wider(
      names_from  = Method,
      values_from = MAE
    )
  
  # MAE1 pair
  pairs_MAE1 <- mae_wide %>%
    filter(!is.na(MAE1a), !is.na(MAE1b)) %>%
    transmute(
      Dataset,
      Sample,
      Marker,
      pair        = "MAE1",
      MAE_emp     = MAE1a,
      MAE_bayes   = MAE1b
    )
  
  # MAE2 pair
  pairs_MAE2 <- mae_wide %>%
    filter(!is.na(MAE2a), !is.na(MAE2b)) %>%
    transmute(
      Dataset,
      Sample,
      Marker,
      pair        = "MAE2",
      MAE_emp     = MAE2a,
      MAE_bayes   = MAE2b
    )
  
  # combine
  pairs_all <- bind_rows(pairs_MAE1, pairs_MAE2)
  
  # make pair an ordered factor to control facet column order
  pairs_all <- pairs_all %>%mutate(pair = factor(pair, levels = c("MAE1", "MAE2")))
  pairs_all$Dataset = factor(pairs_all$Dataset, levels=datNames)
  
  maxLim = 4
  p_shrink <- ggplot(pairs_all, aes(x = MAE_emp, y = MAE_bayes)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.4, linetype = "dashed") +  # regression fit
    facet_grid(Dataset ~ pair) +
    labs(
      x = "Empirical MAE estimate (a)",
      y = "Bayesian MAE estimate (b)",
      title = "Shrinkage patterns of MAE estimators across datasets"
    ) + theme_bw() +  
    coord_cartesian(xlim = c(0, maxLim), ylim = c(0, maxLim))   # <-- added here
  
  ggsave(paste0(outfold,"//Figure_ShrinkageEffect.svg"),plot = p_shrink, width = 6,height=6)
}


compareTimes = function(df, outfold) {
  library(data.table)
  setDT(df)
  
  summary_dt <- df[, .(median_time = ceiling(median(Time, na.rm = TRUE))),  by = .(Dataset, Method)]
  write.table(summary_dt,file=paste0(outfold,"//medianTime_all.csv"),sep=";",row.names = FALSE)
  
  # Get median Time per Dataset, Method, and isComb
  summary_dt2 <- df[, .(median_time = ceiling(median(Time, na.rm = TRUE))), by = .(Dataset, Method, isComb)]
  write.table(summary_dt2,file=paste0(outfold,"//medianTime_CombGrp.csv"),sep=";",row.names = FALSE)
}
