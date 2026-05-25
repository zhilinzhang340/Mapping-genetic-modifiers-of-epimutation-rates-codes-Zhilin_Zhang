# ==============================================================================
# Script Name: plot_Figure4D.R
# Description: Generates Figure 4D - expression VIM2/4 vs global rate RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#fig4D VIM4 expression correlate with rate in RIL


rm(list=ls())
combined_df2 <- fread("/mnt/int/RIL/for_paper/codes_github/data/expression_VIM4_vs_global_rate_2026-01.txt")

combined_df2$peak_genotype <- factor(combined_df2$peak_genotype,levels=c("peakLer","peakCvi"))



output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig4D Scatter plot of VIM4_expression vs Average_rate add color of peak genotype 2025-04"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 绘制 VIM2/4 和 Alpha 的散点图并添加回归线
gg <- ggplot(combined_df2, aes(x = VIM4, y = avg)) +
  geom_point(aes(color = peak_genotype)) +  # 根据分组着色点
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加回归线和置信区间
  labs(title = "",
       x = "Expression level of VIM2/4",
       y = "Avearage_rate") +
  theme_classic()+
  stat_regline_equation(label.x = 9.4, label.y = 6) +  # 显示回归方程
  stat_cor(label.x = 9.4, label.y =5, method = "pearson")+  # 显示相关系数
  scale_color_manual(values = c("peakCvi" = "#00A29A", "peakLer" = "#1D2088")) +  # 自定义颜色
  theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5))

print(gg)
dev.off()









