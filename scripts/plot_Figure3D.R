# ==============================================================================
# Script Name: plot_Figure3D.R
# Description: Generates Figure 3D - dot plot alpha rate  Lerseg_withoutchr1 vs Cviseg_withoutchr1 RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#new RIL main figures fig3d gain loss rate cviseg lerseg


rm(list=ls())

combined_dfalls3 <- fread("/mnt/int/RIL/for_paper/codes_github/data/global_gain_loss_rate_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_2026-01.txt")


combined_dfalls3 $condition <- factor(combined_dfalls3 $condition ,levels=c("All","Lerseg","Cviseg"))

combined_dfalls3 $peak_genotype <- factor(combined_dfalls3 $peak_genotype ,levels=c("peakLer","peakCvi"))



output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 6
out.name<-"Fig3D alpha rate dot plot  Lerseg_withoutchr1 vs Cviseg_withoutchr1 RIL merge peakLer peakCvi"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(combined_dfalls3, aes(x = condition, y = Alpha, fill = peak_genotype)) +
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
    label = "p.format",  
    label.y = 5.2,
    size = 4
  )+
  labs(title = "Cvi vs Ler global withoutchr1", x = "", y = expression(paste('Gain Rate, ',alpha,' (x',10^-4,')'))) +
  scale_fill_manual(values = c( "#1D2088", "#00A29A")) +
  theme_classic() +

  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  )

print(gg)
dev.off()



combined_dfalls3$Beta <- combined_dfalls3$Beta/1e-04

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 6
out.name<-"Fig3D Beta rate dot plot Lerseg_withoutchr1 vs Cviseg_withoutchr1 RIL merge peakLer peakCvi"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(combined_dfalls3, aes(x = condition, y = Beta, fill = peak_genotype)) +
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
    label = "p.format",  
    label.y = 22.2,
    size = 4
  )+
  labs(title = "Cvi vs Ler global withoutchr1", x = "", y =expression(paste('Loss Rate, ',beta,' (x',10^-4,')'))) +
  scale_fill_manual(values = c( "#1D2088", "#00A29A")) +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  )

print(gg)
dev.off()
