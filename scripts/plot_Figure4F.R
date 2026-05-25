
# ==============================================================================
# Script Name: plot_Figure4F.R
# Description: Generates Figure 4F - qPCR of VIM2 VIM4
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

rm(list=ls())

library(ggplot2)


df <- data.frame(
  Sample = c("Col-0",  "RIL50_G2_1.1", "RIL50_G2_2.1",  "RIL149_G2_1.1", "RIL149_G2_2.1"),
  FoldChange = c(1.00, 10.80, 10.29,  0.57, 1.18)
)


df$type <- c("#585959","red","red","purple","purple")


output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name<-"Fig4F qPCR of VIM2 VIM4"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

gg <- ggplot(df, aes(x = Sample, y = FoldChange, fill = Sample)) +
  geom_col(show.legend = FALSE, width = 0.7, fill = df$type,alpha=0.8) +
  geom_text(aes(label = round(FoldChange, 1)), 
            vjust = -0.5, size = 3) +
  labs(
    title = "Relative Expression of notVIM3 (qPCR)",
    y = expression(paste("Relative expression of VIM2 & VIM4 by qPCR")),
    x = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


print(gg)
dev.off()
