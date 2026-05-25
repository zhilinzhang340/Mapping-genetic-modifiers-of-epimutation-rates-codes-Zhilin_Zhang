
# ==============================================================================
# Script Name: plot_Figure2H.R
# Description: Generates Figure 2H - gbM methy level cviseg lerseg
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


#new RIL main figure fig2G methy level cviseg lerseg
rm(list=ls())
combined_dfalls3 <- fread("/mnt/int/RIL/for_paper/codes_github/data/gbM_methy_lvl_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_2026-01.txt")

combined_dfalls3 $condition <- factor(combined_dfalls3 $condition ,levels=c("All","Lerseg","Cviseg"))

combined_dfalls3 $peak_genotype <- factor(combined_dfalls3 $peak_genotype ,levels=c("peakLer","peakCvi"))

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig2H gbM dot plot new methy lvl Lerseg_withoutchr1 vs Cviseg_withoutchr1 RIL merge peakLer peakCvi add average parent"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(combined_dfalls3, aes(x = condition, y = mean_methy_lvl, fill = peak_genotype)) +
  stat_summary(
    fun.data = mean_se, geom = "errorbar", width = 0.3,
    position = position_dodge(width = 0.6), color = "black"
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 21, size = 4,
    position = position_dodge(width = 0.6)
  ) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = after_stat(round(y, 1))),
    vjust = -1.2, color = "red",
    position = position_dodge(width = 0.6)
  ) +
  stat_compare_means(
    aes(group = peak_genotype),
    method = "t.test",
    label.y = 34,
    label = "p.format",
    size = 4
  )+
  labs(
    title = "Cvi vs Ler gbM withoutchr1",
    x = "",
    y = "new methy lvl (%)"
  ) +
  scale_fill_manual(values = c( "#1D2088", "#00A29A")) +
  theme_classic() +
  #coord_cartesian(ylim = c(6, 32)) +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  )

print(gg)
dev.off()





