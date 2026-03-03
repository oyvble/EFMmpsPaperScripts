#script to compare degradation estimates obtained from different methods (compared against default)
#rm(list=ls())
folds  = c("Results","Supplementary")
methods0 = c("MAE1a","MAE1b","MAE2a","MAE2b")  #methods to consider
#methods1 = paste0(methods0,"_ND")  #methods to consider
#methods = c(methods0,methods1)

df = NULL
for(fold in folds) {
  df1 = readRDS(file=paste0(fold,"/df_SummaryResults.RDS"))
  df = rbind(df,df1)
}


metrics = c("Accuracy","Brier","Acc_high","Coverage","Markers95")
library(dplyr)
library(tidyr)

maxDiffTable = NULL
for(met in methods0) {
#  met = methods0[1]
  met2 = paste0(met,"_ND")

  # 1) Diff per Group x Metric (Option B)
  cmp_by_group <- df %>%
    filter(Method %in% c(met, met2), isCombined == "All") %>%
    select(Method, Group, all_of(metrics)) %>%
    pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
    group_by(Group, Metric) %>%
    summarise(
      val_met  = Value[Method == met][1],
      val_met2 = Value[Method == met2][1],
      Diff     = val_met - val_met2,
      .groups = "drop"
    )
  
  # 2) Aggregate across Group -> one value per Metric using max |Diff|
  cmp_agg <- cmp_by_group %>%
    group_by(Metric) %>%
    slice_max(order_by = abs(Diff), n = 1, with_ties = FALSE) %>%  # keep the group that drives max |diff|
    transmute(
      Metric,
      max_abs_Diff = abs(Diff),
      Diff_at_max  = Diff,
      Group_at_max = Group
    ) %>%
    arrange(match(Metric, metrics))
  
  val = setNames(cmp_agg$max_abs_Diff,cmp_agg$Metric)
  maxDiffTable = rbind(maxDiffTable,val)
}
rownames(maxDiffTable) = methods0
outtab = round(maxDiffTable,5)

#outFold
write.table(outtab,file=paste0(outfold,"//comparedMetrics_noDegApproach.csv"),sep=";")
