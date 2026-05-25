
# ==============================================================================
# Script Name: plot_Figure5C.R
# Description: Generates Figure 5C - natural_accessions_global_Methylation_level_across_SNPs_within_deletion
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================
#global_methy_lvl
rm(list=ls())

library(data.table)
library(ggplot2)
library(ggpubr)

usec <- fread("/mnt/int/RIL/for_paper/Table/TableSNew_SNPs_within_deletion_global_mCG_level_across_799natural_accessions.txt")


# 确保数据是 data.table 格式
usec <- as.data.table(usec[,-1])

# 获取 SNP 列名（从第2列到最后）
snp_cols <- names(usec)[2:ncol(usec)]

# 初始化结果列表
result_list <- list()

# 遍历每个 SNP 列
for (snp in snp_cols) {
  #snp <- "S24586079"
  # 创建一个临时表，只含该 SNP 和 global_methy_lvl
  tmp <- usec[, .(SNP = snp, 
                  group = get(snp), 
                  global_methy_lvl)]
  
  # 按 major / minor 分组求平均
  stat <- tmp[, .(mean_methy = mean(global_methy_lvl, na.rm = TRUE)), by = .(SNP, group)]
  
  # 转宽格式：一行一个 SNP，列为 major 和 minor 的均值
  stat_wide <- dcast(stat, SNP ~ group, value.var = "mean_methy")
  
  # 添加到列表
  result_list[[snp]] <- stat_wide
}



# 合并为一个 data.table
snp_summary <- rbindlist(result_list, fill = TRUE)


snp_summary <- snp_summary[2:13,]



# 假设 snp_summary 是 data.table，如果不是先转换
snp_summary <- as.data.table(snp_summary)

# 转成长格式
snp_long <- melt(snp_summary, 
                 id.vars = "SNP", 
                 variable.name = "Group", 
                 value.name = "Methylation")

# 把 Group 设为 factor，顺序为 major / minor
snp_long[, Group := factor(Group, levels = c("major", "minor"))]



snp_long
snp_long$Methylation <- 100*snp_long$Methylation


output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name <-"Fig5C global Methylation level across SNPs within deletion by major minor only SNP within deletion"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 画图
# 绘制 boxplot + jitter（散点）
gg <-ggplot(snp_long, aes(x = Group, y = Methylation, color = Group)) +
  geom_boxplot(width = 0.5,outlier.shape = NA,  fill = NA) +  # 只画箱线不填充
  geom_jitter(width = 0.2,size = 2,  alpha = 0.8) +            # 添加散点
  scale_color_manual(values = c("major" = "#4DBBD5", "minor" = "#E64B35")) +
  stat_compare_means(method = "t.test", label.y = 26.5) +  # 添加 t 检验 p 值
  labs(title = "natural accessions",
       x = "Genotype group",
       y = "Mean global methylation level (per SNP)") +
  theme_classic() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")


print(gg)
dev.off()

