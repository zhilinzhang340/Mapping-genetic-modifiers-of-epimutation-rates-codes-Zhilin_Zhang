
# ==============================================================================
# Script Name: plot_Figure4G.R
# Description: Generates Figure 4G - heatmap correlation gene expression
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

rm(list=ls())
library("pheatmap")

combined_df2 <- fread("/mnt/int/RIL/for_paper/codes_github/data/preTableSNew_expression_level_across_67RILs.txt")



combined_df20 <- combined_df2[combined_df2$peak_genotype=="peakCvi",]
combined_df22 <- combined_df2[combined_df2$peak_genotype=="peakLer",]


usesler <- combined_df22[,-c(1,21)]
head(usesler)
# 计算相关性矩阵
cor_matrix <- cor(usesler)



col_names_ordered <-row_names_ordered <- c("Vim4","MET1","Vim1","CMT3","Vim3","Vim5","Vim6","CMT1","CMT2","DRM1","DRM2",
                                           "DML2","DML3","DME","ROS1","ROS3","ROS4","ROS5","REF6")




# 计算 usescvi 的相关性矩阵
cor_matrix2 <- cor(usesler)

# 根据前一个 heatmap 的顺序对行列重新排序
cor_matrix2_ordered_ler <- cor_matrix2[row_names_ordered, col_names_ordered]


# 
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 8
plot.width <- 13
out.name <- "Fig4G_only peakler_heatmap correlation VIM gene family and MET1 expression and other important genes in 46 RILlines add more genes reorder"

# 定义颜色渐变范围，使 0 为白色，负数为蓝色，正数为红色，1 为红色的最高值
breaks <- seq(-1, 1, length.out = 50) # 根据 -1 到 1 设置分段


p <-pheatmap(cor_matrix2_ordered_ler,
             display_numbers = TRUE, # 在热图上显示相关系数数值
             color = colorRampPalette(c("blue", "white", "red"))(50), # 蓝-白-红的颜色渐变
             breaks = breaks, # 使用自定义的颜色分段
             #cluster_rows = FALSE, # 不进行行聚类
             #cluster_cols = FALSE, # 不进行列聚类
             main = "Correlation of Genes expression in 46 RILlines peakler",
             filename = paste(output.dir, out.name, ".pdf", sep = ""), # 直接保存为 PDF
             width = plot.width,
             cluster_rows = F,
             cluster_cols =F,
             height = plot.height)
print(p)
dev.off()



usescvi <- combined_df20[,-c(1,21)]
head(usescvi)
# 计算相关性矩阵
cor_matrix <- cor(usescvi)



#reorder keep 67RIL order

# 计算 usescvi 的相关性矩阵
cor_matrix2 <- cor(usescvi)

# 根据前一个 heatmap 的顺序对行列重新排序
cor_matrix2_ordered_cvi <- cor_matrix2[row_names_ordered, col_names_ordered]


# 
# 设置输出目录和文件名
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 8
plot.width <- 13
out.name <- "Fig4G_only peakCvi heatmap correlation VIM gene family and MET1 expression and other important genes in 21 RILlines add more genes reorder"

# 定义颜色渐变范围，使 0 为白色，负数为蓝色，正数为红色，1 为红色的最高值
breaks <- seq(-1, 1, length.out = 50) # 根据 -1 到 1 设置分段

# 使用 pheatmap 的 filename 参数直接保存 PDF
pheatmap(cor_matrix2_ordered_cvi ,
         display_numbers = TRUE, # 在热图上显示相关系数数值
         color = colorRampPalette(c("blue", "white", "red"))(50), # 蓝-白-红的颜色渐变
         breaks = breaks, # 使用自定义的颜色分段
         #cluster_rows = FALSE, # 不进行行聚类
         #cluster_cols = FALSE, # 不进行列聚类
         main = "Correlation of Genes expression in 21 RILlines peakCvi",
         filename = paste(output.dir, out.name, ".pdf", sep = ""), # 直接保存为 PDF
         width = plot.width,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         height = plot.height)
