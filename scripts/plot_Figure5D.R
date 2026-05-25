

# ==============================================================================
# Script Name: plot_Figure5D.R
# Description: Generates Figure 5D - natural_accessions_global_Methylation_level_variance_across_SNPs_within_deletion
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================
#global_methy_lvl var

# 清空环境并设置目录
rm(list=ls())

cout <- fread("/mnt/int/RIL/for_paper/codes_github/data/natural_accessions_global_Methylation_level_variance_across_SNPs_within_deletion_by_major_minor.txt")

cout <- cout[3:26,]
cout

output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name <-"Fig5D global Methylation level variance across SNPs within deletion by major minor only SNP within deletion"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 画图
# 绘制 boxplot + jitter（散点）
gg <- ggplot(cout, aes(x = genotype_column, y = variance, color = genotype_column)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA) +  # 只画箱线不填充
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +            # 添加散点
  scale_color_manual(values = c("major" = "#4DBBD5", "minor" = "#E64B35")) +
  labs(title = "natural accessions",
       x = "Genotype group",
       y = "Variance of global methylation level (per SNP)") +
  theme_classic() +
  stat_compare_means(method = "t.test", label.y = 3.2) +  # 添加 t 检验 p 值
  scale_y_continuous(limits = c(1.8,3.4))+
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")


print(gg)
dev.off()

