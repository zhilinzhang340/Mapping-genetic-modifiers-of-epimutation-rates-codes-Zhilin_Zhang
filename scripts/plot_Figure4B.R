
# ==============================================================================
# Script Name: plot_Figure4B.R
# Description: Generates Figure 4B - TE methy AT1TE80755 vs expression VIM2/4 RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


#check new TE methy AT1TE80755 vs expression RIL


rm(list=ls())
combined_df2 <- fread("/mnt/int/RIL/for_paper/codes_github/data/AT1TE80755_methy_lvl_vs_expression_VIM4_2026-01.txt")


combined_df2$peak_genotype <- factor(combined_df2$peak_genotype,levels=c("peakLer","peakCvi"))

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig4B Scatter plot of TE: AT1TE80755 methy_lvl vs VIM4  with regression line RIL add color of peak genotype"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)



gg <- ggplot(combined_df2, aes(x = mean_meth_lvl, y = Vim4)) +
  geom_point(aes(color = peak_genotype)) +  # 根据分组着色点
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # 一条统一的回归线
  labs(
    y = "Expression level of VIM2/4",
    x = "Methylation level of TE: AT1TE80755") +
  theme_classic() +
  stat_regline_equation(label.x = 0.22, label.y = 9.67) +  # 显示回归方程
  stat_cor(label.x = 0.22, label.y = 9.6, method = "pearson") +  # 显示相关系数
  scale_color_manual(values = c("peakCvi" = "#00A29A", "peakLer" = "#1D2088")) +  # 自定义颜色
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))


print(gg)
dev.off()
