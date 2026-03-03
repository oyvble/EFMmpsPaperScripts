  

#Function that creates results for a study
createStudyResults = function(df, dfSummary, ruleList, outFold="Results", highConfidence = 0.95) {
  library(data.table)
  library(flextable);library(officer)
  setDT(df)  
  
  #Helpfunction to save a table as docx
  saveTable = function(table,fn="",header="") {
    ft <- flextable(table) |>
      add_header_lines(values = header) |> theme_booktabs() |> autofit()
    doc <- read_docx() |>  body_add_flextable(ft)
    print(doc, target = paste0(outFold,"/",fn,".docx"))
  }
  
  #Assign rules to the dfs:
  rules = c("Easy","Medium","Hard")
  df$Group=NULL
  dfSummary$Group=NULL
  for(evid in names(ruleList)) {
  #  evid = names(ruleList)[1]
    ruleListEvid = ruleList[[evid]]
    for(ref in names(ruleListEvid)) {
# ref = names(ruleListEvid)[1]
      rule = ruleListEvid[[ref]]
      df$Group[df$Sample==evid & df$Ref==ref] = rule
      dfSummary$Group[dfSummary$Sample==evid & dfSummary$Ref==ref] = rule
    }
  }
  df$Group = factor(df$Group,rules) #make as factor
  dfSummary$Group = factor(dfSummary$Group,rules) #make as factor
  dfSummary$isCombined = grepl("/",dfSummary$Sample)
#  View(dfSummary)
  
  #Obtain summary table (aggregated number of markers and predicted)
  setDT(dfSummary)
  summary_long <- dfSummary[, .(
    n              = .N,
    avg = round(mean(nMarkers, na.rm = TRUE)),
    #total = sum(nMarkers, na.rm = TRUE),
    #total_nMatches = sum(nMatches, na.rm = TRUE),
    mean_accuracy  = round(sum(nMatches, na.rm = TRUE) / sum(nMarkers, na.rm = TRUE)*100)
  ), by = .(isCombined, Group, Method)]
  summary_first <- summary_long[Method == "default",]
  summary_first[, c("mean_accuracy", "Method") := NULL]
  acc_wide <- dcast(summary_long, isCombined + Group ~ Method,value.var = "mean_accuracy")
  final <- merge(summary_first, acc_wide, by = c("isCombined","Group"), all.x = TRUE)[order(isCombined, Group)]
  saveTable(final,fn="Table_Groupoverview",header="Group table") #store
  
  
  #Obtain summary table (obtain best accuracy per method)
  setDT(dfSummary) 
  dt_best <- dfSummary[, .SD[which.max(Acc)], by = .(Sample, Ref)]
  
  summaryTable = dt_best[,.(Sample,Ref,Mx,Drop,Group,Acc)]
  summaryTable$Acc = round(summaryTable$Acc,2)
  summaryTable$Mx = round(summaryTable$Mx,2)
  summaryTable$Drop = round(summaryTable$Drop,2)
  saveTable(summaryTable,fn="Table_Sampleoverview",header="Sample table") #store
  
  #Show in boxplot
  library(ggplot2)
  df_acc = summaryTable
  df_acc$isComb = grepl("/",df_acc$Sample)
  gg = ggplot(df_acc, aes(x = Group, y = Acc, fill = isComb)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_jitter(color = "black", position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.6, size = 1.5) +
    labs( x = "Group", y = "Accuracy", fill = "Combined?", title = "Top Accuracy (among all methods)")
  ggsave(filename = paste0(outFold,"/Figure_TopAccWithGroup.svg"),plot=gg,width=6,height=6)
  
    

  #threshold for high-confidence calls
  #df0 = df; df=df0
  df$isCombined = grepl("/",df$Sample)
  hasCombined = any(df$isCombined)
  
  #RENAMING METHOD NAMES (POSSIBLE)
  if(0) {
    rename0 = c("default","binned") #unchanged
    rename1 = c("MAE1a","MAE1b","MAE2a","MAE2b") #unchanged
    renaming = c(setNames(rename0,rename0),setNames(rename1,rename1))
    #renaming = c(renaming, setNames(c("MAE1a","MAE1b","MAE2a","MAE2b"),c("MAEnaive","MAE","MAE2naive","MAE2")))
    df$Method = renaming[df$Method]
    #table(df$Method)
  }  
  # 1) Keep ONLY the top-ranked genotype per Method×Sample×Ref×Locus
  #Get top probs, then take the mean match of these = accuracy
  top <- df[order(Method, Sample, Ref, Locus, -Prob, -Match)][, .SD[1L], by = .(Method, Sample, Ref, Locus)]
  acc_by_m <- top[, .(Acc = mean(Match)), by = .(Sample, Ref, Method,Group,isCombined)]
  
  #ACCURACY
  overallAcc <- acc_by_m[, .(Accuracy = mean(Acc)), by = .(Method, Group)]
  overallAcc[, isCombined := "All"]  # label to align with split table
  if(hasCombined) {
    overallAcc_Combined <- acc_by_m[, .(Accuracy = mean(Acc)), by = .(Method, Group,isCombined)]
    overallAcc <- rbindlist(list(overallAcc, overallAcc_Combined), use.names = TRUE, fill = TRUE)
  }
  
  #BRIER SCORE
  overallBrier <- df[, .(Brier = mean((Prob - Match)^2)), by = .(Method, Group)]
  overallBrier[, isCombined := "All"]  # label to align with split table
  if(hasCombined) {
    overallBrier_Combined <- df[, .(Brier = mean((Prob - Match)^2)), by = .(Method, Group, isCombined)]
    overallBrier <- rbindlist(list(overallBrier, overallBrier_Combined), use.names = TRUE, fill = TRUE)
  }
  overall <- overallAcc[overallBrier, on = .(Method, Group,isCombined)] # merge tables
  
  
  #4) Restrict predictions with p>0.95. Obtain accuracy/coverage/Brier
  
  # .0) collapse to ONE row per marker = top-ranked genotype
  #    tie-breaker: prefer higher Prob; if tie, prefer correct (Match=1)
  top <- df[order(Method, Group, isCombined, Sample, Ref, Locus, -Prob, -Match)
  ][, .SD[1L], by = .(Method, Group, isCombined, Sample, Ref, Locus)]
  
  # .1) Totals = number of markers per group (denominator for coverage)
  totals <- top[, .(TotalMarkers = .N), by = .(Method, Group, isCombined)]
  #sum(totals$TotalMarkers)/134
  
  # .2) High-confidence subset: markers where top-call Prob >= 0.95
  high <- top[Prob >= highConfidence,
              .(Acc_high   = mean(Match), Brier_high = mean((Prob - Match)^2), Markers95  = .N),
              by = .(Method, Group, isCombined)]
  
  # 3) Merge & compute coverage
  summary95 <- high[totals, on = .(Method, Group, isCombined)]
  summary95[, Coverage := Markers95 / TotalMarkers]
  
  if(hasCombined) {
    #Also compute “overall” (ignoring isCombined) from the same top table
    totals_all <- top[, .(TotalMarkers = .N), by = .(Method, Group)]
    high_all   <- top[Prob >= highConfidence,
                      .(Acc_high   = mean(Match), Brier_high = mean((Prob - Match)^2), Markers95  = .N),
                      by = .(Method, Group)]
    summary95_all <- high_all[totals_all, on = .(Method, Group)]
    summary95_all[, Coverage := Markers95 / TotalMarkers]
    summary95_all[, isCombined := "All"]  # label to align with split table
    summary95 <- rbindlist(list(summary95_all, summary95), use.names = TRUE, fill = TRUE)     #Combine overall + split-by-isCombined
    final = overall[summary95, on = .(Method, Group, isCombined)] #merge overall + high-confidence results
  } else {
    final = overall[summary95, on = .(Method, Group)] #merge overall + high-confidence results
    final$i.isCombined = NULL
  }
  #  View(final)
  saveRDS(final,file=paste0(outFold,"/df_SummaryResults.RDS"))
  
  final$isCombined = factor(final$isCombined,levels=c("TRUE","FALSE","All"))
  modifyZero = function(x) gsub("0\\.","\\.",sprintf("%.2f", x))
  #modifyZero(final$Brier)
  #Round numeric columns and order nicely
  colPerc = grepl("Acc|Coverage",names(final))
  final[, (names(final)[colPerc]) := lapply(.SD, function(x) paste0(round(x*100))), .SDcols = colPerc]
  num_cols <- grepl("Brier",names(final))
  final[, (names(final)[num_cols]) := lapply(.SD, modifyZero), .SDcols = num_cols]
  final = final[order(isCombined,Group, Method)][, .SD[1L], by = .(Group, isCombined, Method)]
  final = final[, .(isCombined,Group, Method, Brier, Accuracy, Acc_high, Coverage, Markers95)]#,Calibrated)]
#  View(final)
  #Brier_high, Markers95, TotalMarkers,
  #if(length(unique(final$isCombined))>1 && removeAll) final = final[final$isCombined!="All",]
  #Create summary table
  saveTable(final,fn="Table_SummaryResults",header="Comparison")
  
  
  #4) Showing calibration (full range):
  #sdf = df[!df$isCombined,]
  makeCalibrationPlot = function(sdf,suffix="",breaks=NULL) {
    #sdf = df  
    if(is.null(breaks)) breaks <- seq(0, 1, length.out = 21)
    #breaks <- seq(0, 1, length.out = 11)
    sdf[, Bin := cut(Prob, breaks = breaks, include.lowest = TRUE)]
    
    # 4) Bin-level summaries per Method×Group
    by_bin <- sdf[, .(
      n      = .N,
      mean_p = mean(Prob), # mean predicted probability in bin
      obs    = mean(Match) # observed frequency correct in bin
    ), by = .(Method, Group, Bin)]
    by_bin <- by_bin[!is.na(mean_p) & !is.na(obs)]
    # Note: by_bin must have: Method, Group, Bin, n, mean_p, obs
  
    # Adding 95% Wilson confidence interval for observed frequency per bin ---
    # ChatGPT proposed this: Wilson interval is good for small n
    wilson_ci <- function(k, n, conf = 0.95){
      if(n == 0) return(c(NA, NA))
      z <- qnorm(1 - (1 - conf)/2)
      p <- k/n
      denom <- 1 + z^2/n
      centre <- (p + z^2/(2*n)) / denom
      half <- (z * sqrt( (p*(1-p)/n) + (z^2/(4*n^2)) )) / denom
      c(centre - half, centre + half)
    }
    
    bb <- copy(by_bin)
    bb[, `:=`(
      k = round(obs * n),  # successes
      ci_low = NA_real_, ci_high = NA_real_
    )]
    
    bb[, c("ci_low","ci_high") := {
      ci <- wilson_ci(k, n)
      list(ci[1], ci[2])
    }, by = .(Method, Group, Bin)]
    
    # --- Calibration plot: obs vs mean_p, with diagonal y=x ---
    gg =  ggplot(bb, aes(x = mean_p, y = obs)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      #geom_abline(slope = 1, intercept = 0.01, linetype = 3) +
      #geom_abline(slope = 1, intercept = -0.01, linetype = 3) +
      geom_abline(slope = 1, intercept = 0.05, linetype = 2,col=2) +
      geom_abline(slope = 1, intercept = -0.05, linetype = 2,col=2) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.0, alpha = 0.5) +
      geom_point(alpha = 0.8) +
      # Optional smooth (visual guide only):
      geom_smooth(method = "loess", se = FALSE, formula = y ~ x, span = 0.8,linewidth = 0.5) +
      #scale_size_continuous(name = "Bin size (n)") +
      coord_equal(xlim = c(0,1), ylim = c(0,1), expand = TRUE) +
      labs(title = "Reliability (Calibration) by Method and Group",
           x = "Mean predicted probability in bin",
           y = "Observed frequency correct in bin") +
      facet_grid(Group~Method) + theme(legend.position = "none")
    
    ggsave(filename = paste0(outFold,"/Figure_CalibrationPlot",suffix,".svg"),plot=gg,width=12,height=6)
  }
  makeCalibrationPlot(df,suffix="_all")
  if(hasCombined) {
    makeCalibrationPlot(df[!df$isCombined,],suffix="_NoComb")
    makeCalibrationPlot(df[df$isCombined,],suffix="_Comb")
  }
}



