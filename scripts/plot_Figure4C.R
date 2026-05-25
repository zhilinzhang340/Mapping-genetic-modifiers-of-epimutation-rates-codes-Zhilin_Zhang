

# ==============================================================================
# Script Name: plot_Figure4C.R
# Description: Generates Figure 4C - expression VIM2/4 vs global methylation level RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#new fig4C VIM2/4 expression with methylation level in RIL data 2025-05


rm(list=ls())
combined_df2 <- fread("/mnt/int/RIL/for_paper/codes_github/data/expression_VIM4_vs_global_methy_lvl_2026-01.txt")

combined_df2$peak_genotype <- factor(combined_df2$peak_genotype,levels=c("peakLer","peakCvi"))

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig4C Scatter plot of VIM4_expression vs global methy add color of peak genotype 2025-04"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 绘制 VIM2/4 和 Alpha 的散点图并添加回归线
gg <- ggplot(combined_df2, aes(x = VIM4, y = global_methy_lvl)) +
  geom_point(aes(color = peak_genotype)) +  # 根据分组着色点
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # 添加回归线和置信区间
  labs(title = "",
       x = "Expression level of VIM2/4",
       y = "global methy lvl") +
  theme_classic()+
  stat_regline_equation(label.x = 9.1, label.y = 16) +  # 显示回归方程
  stat_cor(label.x = 9.1, label.y = 15.7, method = "pearson")  +# 显示相关系数
  scale_color_manual(values = c("peakCvi" = "#00A29A", "peakLer" = "#1D2088")) +  # 自定义颜色
  theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5))

print(gg)
dev.off()








