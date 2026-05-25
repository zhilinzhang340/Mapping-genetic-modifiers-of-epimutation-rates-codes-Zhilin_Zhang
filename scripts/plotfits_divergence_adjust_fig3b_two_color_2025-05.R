test_plotFITS <- function(pedigree.names, input.dir, output.dir, out.name, 
                          alpha, geom.point.size, geom.line.size, plot.height, plot.width, plot.type, lsq.line, facets.ncols) {
  
  library(ggplot2)
  
  df_div <- c()
  theory.df_fits <- c()
  all.miss <- c()
  gg <- ggplot()
  
  # 十个 RIL line：颜色 × alpha 组合唯一
  ril_names <- c("RIL21", "RIL58", "RIL149", "RIL150", "RIL182", "RIL50", "RIL15", "RIL155", "RIL162", "RIL7")
  
  ril_colors <- c(
    "RIL21" ="purple", "RIL58" ="purple", "RIL149" ="purple", "RIL150" ="purple", "RIL182" ="purple",
    "RIL50" = "red", "RIL15" = "red", "RIL155" = "red", "RIL162" ="red", "RIL7" ="red"
  )
  
  ril_alpha <- c(
    "RIL21" = 1.0, "RIL58" = 0.8, "RIL149" = 0.6, "RIL150" = 0.5, "RIL182" = 0.4,
    "RIL50" = 1.0, "RIL15" = 0.8, "RIL155" = 0.6, "RIL162" = 0.5, "RIL7" = 0.4
  )
  
  for (i in 1:NROW(pedigree.names)) {
    for (j in 1:NCOL(pedigree.names)) {
      
      datain <- dget(paste0(input.dir, pedigree.names[i, j]))
      
      name <- gsub("(_).*_ABneutral_CG_estimates.Rdata$", "", basename(as.character(pedigree.names[i, j])))
      Genotypes <- gsub(pattern = paste0(name, "_|\\_ABneutral_CG_estimates.Rdata$"), "", 
                        basename(as.character(pedigree.names[i, j])))
      
      RIL <- name
      intercept <- datain$estimates[1, "intercept"]
      delta.t <- datain$pedigree[, "delta.t"]
      div.obs <- datain$pedigree[, "div.obs"] - intercept
      df_div <- rbind(df_div, data.frame(Genotypes, RIL, delta.t, div.obs))
      
      miss <- cbind(paste0(RIL, "_", Genotypes), length(div.obs[div.obs < 0]))
      all.miss <- rbind(all.miss, miss)
      
      theory.fit.data <- datain$for.fit.plot
      theory.fits <- c(intercept, theory.fit.data[, "div.sim"]) - intercept
      theory.fit.t <- c(0, theory.fit.data[, "delta.t"])
      theory.df_fits <- rbind(theory.df_fits, data.frame(Genotypes, RIL, theory.fit.t, theory.fits))
    }
  }
  
  write.table(all.miss, "output_miss_point_number_below_zero_during_divergence_plot.txt",
              append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ymax <- 1.1 * max(as.numeric(df_div$div.obs), na.rm = TRUE)
  xmax <- 1.01 * max(as.numeric(df_div$delta.t), na.rm = TRUE)
  
  Genotypes_1 <- strsplit(as.character(df_div$Genotypes[1]), "_")[[1]]
  Genotypes_2 <- paste(Genotypes_1, collapse = " ")
  
  pdf(paste0(output.dir, out.name, ".pdf"), width = plot.width, height = plot.height)
  
  gg <- gg +
    geom_point(data = df_div, aes(x = delta.t, y = div.obs, color = RIL),
               alpha = alpha, size = geom.point.size) +
    geom_line(data = theory.df_fits, aes(x = theory.fit.t, y = theory.fits, color = RIL, alpha = RIL),
              size = geom.line.size, linetype = "solid") +
    labs(
      x = expression(paste(Delta, "t (generations)")),
      y = "mCG divergence",
      color = "RIL Line",
      alpha = "RIL Line",
      title = Genotypes_2
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    scale_color_manual(values = ril_colors) +
    scale_alpha_manual(values = ril_alpha) +
    guides(
      color = guide_legend(title = "RIL Line"),
      alpha = guide_legend(title = "RIL Line", override.aes = list(linetype = "solid", size = geom.line.size))
    )
  
  print(gg)
  dev.off()
}
