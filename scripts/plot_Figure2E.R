

#Fig2E
# ==============================================================================
# Script Name: plot_Figure2E.R
# Description: Generates Figure 2E - global_per5kb methylation level peak marker
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#plot
rm(list=ls())
nout <- fread("/mnt/int/RIL/for_paper/codes_github/data/68RIL_global_per5kb_qtl_peak_marker_with_pos.txt")

use <- nout[,c("chr","midpos","gchr","gposmid")]



library(ggplot2)
library(forcats)


# 将 gchr 转换为因子并反转其级别
use$gchr <- fct_rev(as.factor(use$gchr))


output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 10
plot.width <- 10
out.name<-"Fig2E global_per5kb methylation level peak marker by Chromosome"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 创建散点图
gg <-ggplot(use, aes(x = midpos/1e+07, y = gposmid/1e+07)) +
  geom_point() +  # 添加散点
  facet_grid(gchr ~ chr, scales = "free") +  # 根据 chr 和 gchr 创建网格布局
  labs(x = expression(paste("Peak marker Position (x",10^7,")")), y = expression(paste("global marker Position (x",10^7,")")), title = "global_per5kb methylation level by Chromosomes") +
  theme_classic()+
  theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5),)

print(gg)
dev.off()
