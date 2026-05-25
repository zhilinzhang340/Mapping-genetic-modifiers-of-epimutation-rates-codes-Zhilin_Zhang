# ==============================================================================
# Script Name: Global_PCA_and_Target_DEG_Barplot_Separate_PDFs.R
# Description: Normalization, global PCA, limma DEG, and customized target gene barplot (P-value only)
# Author: Zhilin Zhang
# ==============================================================================

rm(list=ls())
gc()

library(data.table)
library(dplyr)
library(ggplot2)
library(limma)
library(stringr)

# 1. 读取原始数据与基因型
# ------------------------------------------------------------------------------
expr_full <- fread("/mnt/int/RIL/for_paper/Table/eQTL_raw_data_67RIL_RIL186miss_all_eqtl_22644_overlap_TAIR10_genes_list.txt")
peak <- fread("/mnt/int/RIL/qtl/eqtl/peak_status_68RIL_1_2623_0as-1_2as1.txt")

# 转换表达矩阵格式
expr_mat <- as.matrix(expr_full[, -1, with=FALSE])
rownames(expr_mat) <- expr_full$RILRIL

# 样本对齐：剔除缺失的 RIL186，确保数据严格匹配
common_rils <- intersect(colnames(expr_mat), peak$RIL)
expr_mat <- expr_mat[, common_rils]
peak_sub <- peak[match(common_rils, peak$RIL), ]

# ==============================================================================
# 基因型映射 (0as-1 -> -1 代表 Cvi, 2as1 -> 1 代表 Ler)
# 设置 levels 确保 Ler 为对照组 (Base level)，计算 LogFC 时为 Cvi - Ler
# ==============================================================================
group <- factor(ifelse(peak_sub$peak_status == 1, "peakLer", "peakCvi"), levels = c("peakLer", "peakCvi"))

# 2. 数据标准化 (Normalization)
# ------------------------------------------------------------------------------
#expr_norm <- normalizeBetweenArrays(expr_mat, method = "quantile")
expr_norm <-expr_mat

# 3. 全转录组 PCA 分析
# ------------------------------------------------------------------------------
pca_res <- prcomp(t(expr_norm), scale. = TRUE)
pca_df <- data.frame(Sample = rownames(pca_res$x), 
                     PC1 = pca_res$x[,1], 
                     PC2 = pca_res$x[,2], 
                     Genotype = group)

pc1_var <- round(summary(pca_res)$importance[2,1] * 100, 1)
pc2_var <- round(summary(pca_res)$importance[2,2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Genotype)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Genotype), geom = "polygon", alpha = 0.1, type = "norm") +
  scale_color_manual(values = c("peakLer" = "#1D2088", "peakCvi" = "#00A29A")) +
  scale_fill_manual(values = c("peakLer" = "#1D2088", "peakCvi" = "#00A29A")) +
  theme_classic() +
  labs(title = "Global PCA (22,644 Genes)", 
       x = paste0("PC1 (", pc1_var, "%)"), 
       y = paste0("PC2 (", pc2_var, "%)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# 4. DEG 分析 (Limma)
# ------------------------------------------------------------------------------
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Cvi_vs_Ler")

fit <- lmFit(expr_norm, design)
fit <- eBayes(fit)
deg_results <- topTable(fit, coef = "Cvi_vs_Ler", number = Inf, adjust.method = "BH")
deg_results$gene_id <- rownames(deg_results)

# 5. 映射并提取靶向表观遗传关键基因
# ------------------------------------------------------------------------------
gene_map <- data.frame(
  gene_id = c("AT5G49160", "AT1G80740", "AT4G19020", "AT1G69770", "AT5G15380", 
              "AT5G14620", "AT1G57820", "AT1G66050", "AT5G39550", "AT1G66040", 
              "AT1G57800", "AT4G08590", "AT5G04560", "AT3G10010", "AT4G34060", 
              "AT2G36490", "AT5G58130", "AT3G14980", "AT1G54840", "AT3G48430"),
  gene_name = c("MET1","CMT1","CMT2","CMT3","DRM1","DRM2",
                "VIM1","VIM2","VIM3","VIM4","VIM5","VIM6",
                "DME","DML2","DML3","ROS1","ROS3","ROS4","ROS5","REF6")
)

res_target <- merge(gene_map, deg_results, by = "gene_id")

# 6. 准备定制化绘图数据 (移除 FDR，仅保留 P.Value)
# ------------------------------------------------------------------------------
plot_data <- res_target %>%
  filter(!is.na(logFC)) %>%
  mutate(
    Type = ifelse(P.Value < 0.05, "Significant (P < 0.05)", "Not Significant")
  )

# 按照截图从上到下的顺序排列
target_order <- c("VIM2", "VIM4", "MET1", "VIM1", "CMT3", "VIM3", "VIM5", "VIM6",
                  "CMT1", "CMT2", "DRM1", "DRM2", "DML2", "DML3", "DME",
                  "ROS1", "ROS3", "ROS4", "ROS5", "REF6")

# 确保提取的数据集中只包含我们要画的基因
plot_data <- plot_data %>% filter(gene_name %in% target_order)

# 反转顺序以适配 coord_flip()
plot_data$gene_name <- factor(plot_data$gene_name, levels = rev(target_order))

# 7. 绘制并分别保存为两个 PDF
# ------------------------------------------------------------------------------
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

# --- 7.1 单独保存 PCA 图 ---
pca.out.name <- "Fig4H_Global_PCA_peakCvi_vs_peakLer"
pdf(paste0(output.dir, pca.out.name, ".pdf"), width = 7, height = 5)
print(p_pca)
dev.off()

# --- 7.2 单独保存 DEG 柱状图 ---
deg_bar <- ggplot(plot_data, aes(x = gene_name, y = logFC, fill = Type)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Significant (P < 0.05)"   = "#CC3333", 
                               "Not Significant"          = "gray70")) +
  labs(
    title = "Differential Expression: peakCvi vs peakLer",
    subtitle = "Log2 Fold Change of Epigenetic Target Genes",
    x = NULL,
    y = "Log2 Fold Change (Positive = Up in peakCvi)",
    fill = "Significance"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 11, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  )

deg.out.name <- "used_Fig4H_Target_Genes_DEG_peakCvi_vs_peakLer"
pdf(paste0(output.dir, deg.out.name, ".pdf"),width = 8, height = 8)
print(deg_bar)
dev.off()

# 保存统计数据
fwrite(plot_data, paste0(output.dir, "used_Fig4H_Target_Genes_DEG_Stats.txt"), sep = "\t")