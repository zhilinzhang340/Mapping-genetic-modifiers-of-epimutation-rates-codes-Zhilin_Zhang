
# ==============================================================================
# Script Name: plot_Figure2I.R
# Description: Generates Figure 2I - global methy level cviseg lerseg only G2 only G9
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================
#new Fig2I dot plot new methy lvl Lerseg_withoutchr1 vs Cviseg_withoutchr1 RIL merge peakLer peakCvi add average parent only G2 only G9
rm(list=ls())
library(ggpubr)
combined_dfalls2 <- fread("/mnt/int/RIL/for_paper/codes_github/data/global_methy_lvl_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_only_G2_only_G9_2026-01.txt")
combined_dfalls2

# 将数字 2 变成 "02"，9 变成 "09"
combined_dfalls2$gen <- sprintf("%02d", combined_dfalls2$gen)


combined_dfalls2$type <- factor(combined_dfalls2$type ,levels = c("Lerseg_peakLer","Lerseg_peakCvi","Cviseg_peakLer","Cviseg_peakCvi")) 



output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig2I dot plot new methy lvl Lerseg_withoutchr1 vs Cviseg_withoutchr1 RIL merge peakLer peakCvi add average parent only G2 only G9"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(combined_dfalls2, aes(x =  type, y = mean_methy_lvl, fill = gen)) +
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
    aes(group = gen),
    method = "t.test",
    label.y = 24,
    label = "p.format",     # 自动保留 2 位有效数字
    size = 4
  ) +
  labs(
    title = "Cvi vs Ler global withoutchr1 only G2 only G9",
    x = "",
    y = "new methy lvl (%)"
  ) +
  #scale_fill_manual(values = c("grey", "purple", "yellow", "cyan", "blue")) +
  theme_classic() +
  scale_fill_manual(values = c("02" = "#F8766D",
                               "09" = "#601986"
  )) +
  scale_color_manual(values = c("02" = "#F8766D",
                                "09" = "#601986"))+
  
  #coord_cartesian(ylim = c(6, 32)) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5)
  )

print(gg)
dev.off()

