
makeSummaryFig = function(df, what,outfold,suffix="all") {
  library(ggplot2)
  
  method_cols <- c(
    "MAE1a"   = "#56B4E9",
    "MAE1b"   = "#009E73",
    "MAE2a"   = "#CC79A7",
    "MAE2b"   = "#D55E00",
    "default" = "#E69F00",
    "binned"  = "#7F7F7F"
  )
  
  #what = "Brier"
  #what = "Coverage"
  yvals = df[[what]]

  df_plot <- df
  df_plot$star_lbl = ifelse(df$Acc_high < 0.94, "*", NA_character_) #indicate if not calibrated
  dodge_w = 0.75
  gg <- ggplot(df_plot, aes(x = Dataset, y = yvals, fill = Method)) +
  geom_col(position = position_dodge(width = dodge_w),width = 0.7,color = "grey30") +
  facet_grid(Group ~ combined, scales = "free", space = "free_x") +
  scale_x_discrete(drop = TRUE) + # make sure discrete x drops unused levels per panel
  scale_fill_manual(values = method_cols, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA)
  )
  if(what=="Brier") gg = gg + labs(y = "Brier score (lower = better)",x = NULL,fill = "Method",title = "Brier score")
  if(what=="Coverage") gg = gg + labs(y = "HC-coverage (higher = better)",x = NULL,fill = "Method",title = "High Confidence Coverage")
  if(what=="Accuracy") gg = gg + labs(y = "Overall Accuracy (higher = better)",x = NULL,fill = "Method",title = "Accuracy for top genotype")
  if(what=="Acc_high") {
    gg = gg + labs(y = "HC-accuracy (higher = better)",x = NULL,fill = "Method",title = "High Confidence Accuracy")
    gg = gg + geom_hline(yintercept = 0.95, color = "black", linetype = "dashed", linewidth = 0.4)
  }
  
  # STAR ANNOTATION: place on top of each bar where HCacc < threshold
  gg <- gg + geom_text( data = df_plot, aes(label = star_lbl), position = position_dodge(width = dodge_w), vjust = +0.4, na.rm = TRUE, size = 4 )
  
  
  # ===== In-panel legend (auto position in top-left facet; top-right corner) =====
  # pick the first facet levels actually present
  if(what=="Brier") {
    g0 <- levels(df_plot$Group)[1]
    c0 <- levels(df_plot$combined)[2]
  } else {
    nGrps = length(unique(df_plot$Group))
    g0 <- levels(df_plot$Group)[nGrps]
    c0 <- levels(df_plot$combined)[2]
  }
  
  y_max = max(yvals[df$Group==g0],na.rm = TRUE) #position legend with this regard
  facet_df <- droplevels(subset(df_plot, Group == g0 & combined == c0))
  if (nrow(facet_df) > 0) {
    # datasets actually present in this facet (respecting dropped levels)
    dsets <- levels(facet_df$Dataset)
    dsets <- dsets[!is.na(dsets)]
    # rightmost dataset index on the discrete scale
    idx_right <- length(dsets)
    
    # compute a safe top y for legend rows
    #y_max <- suppressWarnings(max(yvals, na.rm = TRUE))
    y_min <- suppressWarnings(min(yvals, na.rm = TRUE))
    if (!is.finite(y_min)) y_min <- 0
    if (!is.finite(y_max)) y_max <- 1
    
    # headroom: use existing headroom (from expand) ~10% of range
    yrange <- max(1e-9, y_max - y_min)
    y_top  <- y_max + 0.08 * yrange         # start a bit below very top
    y_step <- 0.08 * yrange                 # vertical spacing between legend rows
    
    # x positions: place boxes/text just to the right of the last dataset's center
    x_box_center <- idx_right #+ 0.35        # tweak for tighter/looser fit
    x_box_w      <- 0.06                    # box half-width ~ discrete units
    x_text       <- x_box_center + 0.10
    
    legend_df <- data.frame(
      Method   = names(method_cols),
      col      = unname(method_cols),
      Group    = factor(g0, levels = levels(df_plot$Group)),
      combined = factor(c0, levels = levels(df_plot$combined)),
      xmin = x_box_center - x_box_w,
      xmax = x_box_center + x_box_w,
      ymin = y_top - (seq_along(method_cols) - 1) * y_step - 0.015 * yrange,
      ymax = y_top - (seq_along(method_cols) - 1) * y_step + 0.015 * yrange,
      xt  = x_text,
      yt  = y_top - (seq_along(method_cols) - 1) * y_step,
      lab = names(method_cols)
    )
    
    gg <- gg +
      ggnewscale::new_scale_fill() +
      geom_rect(
        data = legend_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col),
        inherit.aes = FALSE,
      ) + 
      geom_text(
        data = legend_df,
        aes(x = xt, y = yt, label = lab),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3.1
      )  +  scale_fill_identity(guide = "none") +  # use hex directly, no legend
      coord_cartesian(clip = "off")  
  }
  ggsave(paste0(outfold,"//Figure_allDatasets_",what,"_",suffix,".svg"),plot = gg,width = 10,height=10)
}



