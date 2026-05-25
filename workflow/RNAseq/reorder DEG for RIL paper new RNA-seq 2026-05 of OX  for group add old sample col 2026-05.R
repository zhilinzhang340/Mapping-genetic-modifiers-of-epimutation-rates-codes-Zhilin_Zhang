# ==============================================================================
# 项目：合并所有(26-03, 26-04, 26-05) OX 样本进行 Group DEG 靶基因差异表达分析
# 过滤逻辑：
# 1. 剔除所有 Cvi 样本
# 2. 剔除 OX_Col_VIM4_Col_5-2 样本
# 3. 保留并合并多组纯 Col 背景的 OX 样本 (VIM2, VIM4, VIM24) 与 WT
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 环境初始化与包加载
# ------------------------------------------------------------------------------
rm(list=ls())

library(data.table)
library(tidyverse)
library(DESeq2)

# 设置新的工作目录和输出目录
#setwd("/mnt/ssd/RIL_meth/new_RNA-seq_OX_2026-05/analysis_results/03_counts")
setwd("/mnt/int/RIL/for_paper/codes_github/workflow/RNAseq")

cat("\n=== 1. 开始跨目录批量读取并清理表达矩阵 ===\n")

# ------------------------------------------------------------------------------
# 1. 跨目录批量读取数据、精确过滤样本并构建分组表达矩阵
# ------------------------------------------------------------------------------
# 定义三个目录路径
wt_dir <- "/mnt/ssd/RIL_meth/new_RNA-seq_OX_mutant_2026-03/total_counts"
ox_dir_04 <- "/mnt/ssd/RIL_meth/new_RNA-seq_OX_2026-04/analysis_results/03_counts"
ox_dir_05 <- "/mnt/ssd/RIL_meth/new_RNA-seq_OX_2026-05/analysis_results/03_counts"

# 获取对应目录下的文件
old_files <- list.files(path = wt_dir, pattern = "(Control_Col_WT|OX_).*_counts\\.txt$", full.names = TRUE)
ox_04_files <- list.files(path = ox_dir_04, pattern = "OX_.*_counts\\.txt$", full.names = TRUE)
ox_05_files <- list.files(path = ox_dir_05, pattern = "OX_.*_counts\\.txt$", full.names = TRUE)

# 合并所有文件路径
all_files <- c(old_files, ox_04_files, ox_05_files)

# --- 关键过滤步骤 ---
# 1. 剔除所有文件名中包含 "Cvi" 的样本，保持纯 Col 背景
all_files <- all_files[!grepl("Cvi", basename(all_files), ignore.case = TRUE)]

# 读取文件并构建长表格
long_df <- lapply(all_files, function(f) {
  dt <- fread(f, header = FALSE, col.names = c("gene_id", "count"))
  
  file_name <- basename(f)
  clean_name <- gsub("_counts.txt", "", file_name)
  clean_name <- gsub("OX_Col-", "OX_Col_", clean_name)
  
  dt$sample_name <- clean_name
  return(dt)
}) %>% bind_rows() %>% 
  filter(grepl("^AT[1-5]G", gene_id))

# 转换为宽表格矩阵
counts_df <- long_df %>%
  pivot_wider(names_from = sample_name, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene_id")

counts_matrix <- as.matrix(counts_df)
cat("表达矩阵构建完成！基因数:", nrow(counts_matrix), "，有效样本数:", ncol(counts_matrix), "\n")

# ------------------------------------------------------------------------------
# 2. 构建 Group 实验设计
# ------------------------------------------------------------------------------
cat("\n=== 2. 提取分组名并构建实验设计 ===\n")

sample_names <- colnames(counts_matrix)

# 通过正则去除末尾的重复编号以提取组名 (例如把 _5-1, _7-2 等剥离)
group_names <- sub("_[0-9]+(-[0-9]+)*$", "", sample_names)

colData_group <- data.frame(
  row.names = sample_names,
  Group = factor(group_names)
)

# 强制将野生型设为对比例的 Baseline
colData_group$Group <- relevel(colData_group$Group, ref = "Control_Col_WT")

cat("--- 过滤后的实验分组统计 (包含 2026-05 的新样本) ---\n")
print(table(colData_group$Group))

# ------------------------------------------------------------------------------
# 3. 运行 DESeq2 核心分析
# ------------------------------------------------------------------------------
cat("\n=== 3. 运行 DESeq2 差异分析 ===\n")

dds_group <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                    colData = colData_group,
                                    design = ~ Group)

dds_group <- DESeq(dds_group)

# ------------------------------------------------------------------------------
# 4. 定义 20 个表观遗传核心靶基因
# ------------------------------------------------------------------------------
target_genes <- c("AT5G49160", "AT1G80740", "AT4G19020", "AT1G69770", "AT5G15380", 
                  "AT5G14620", "AT1G57820", "AT1G66050", "AT5G39550", "AT1G66040", 
                  "AT1G57800", "AT4G08590", "AT5G04560", "AT3G10010", "AT4G34060", 
                  "AT2G36490", "AT5G58130", "AT3G14980", "AT1G54840", "AT3G48430")

gene_map <- data.frame(
  gene_id = target_genes,
  gene_name = c("MET1","CMT1","CMT2","CMT3","DRM1","DRM2",
                "VIM1","VIM2","VIM3","VIM4","VIM5","VIM6",
                "DME","DML2","DML3","ROS1","ROS3","ROS4","ROS5","REF6")
)

# ------------------------------------------------------------------------------
# 5. 循环提取每个组别的差异结果并绘制 Barplot
# ------------------------------------------------------------------------------
cat("\n=== 4. 开始提取纯 Col 组靶基因结果并绘制 Barplot ===\n")

# 提取所有的实验组名称 (排除 WT)
treatment_groups <- levels(colData_group$Group)[levels(colData_group$Group) != "Control_Col_WT"]

# 定义固定的靶基因绘图顺序（从上到下）
target_order <- c("VIM2", "VIM4", "MET1", "VIM1", "CMT3", "VIM3", "VIM5", "VIM6",
                  "CMT1", "CMT2", "DRM1", "DRM2", "DML2", "DML3", "DME",
                  "ROS1", "ROS3", "ROS4", "ROS5", "REF6")

for (grp in treatment_groups) {
  
  res <- results(dds_group, contrast = c("Group", grp, "Control_Col_WT"))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # 保存组别的所有基因 DEG 表格
  file_name_all_ox <- paste0("DEGs_All_Genes_", grp, "_vs_WT_NoFilter.txt")
  write.table(res_df, file = file_name_all_ox, quote = FALSE, sep = "\t", row.names = FALSE)
  
  # 匹配关注的 20 个靶基因
  res_target <- merge(gene_map, res_df, by="gene_id")
  
  # 保存靶基因差异表达数据
  txt_file_name <- paste0("DEGs_Target_", grp, "_vs_WT_NoFilter.txt")
  write.table(res_target, file = txt_file_name, quote = FALSE, sep = "\t", row.names = FALSE)
  
  # 整理绘图数据（加入指定的排序逻辑）
  plot_data <- res_target %>%
    filter(!is.na(log2FoldChange)) %>%
    mutate(
      Type = ifelse(!is.na(pvalue) & pvalue < 0.05, "Significant (P < 0.05)", "Not Significant")
    ) %>%
    filter(gene_name %in% target_order)
  
  # 反转顺序以适配 coord_flip() 并在出图时维持与截图相同的上下顺序
  plot_data$gene_name <- factor(plot_data$gene_name, levels = rev(target_order))
  
  plot_title <- paste0(grp, " vs WT")
  plot_y_label <- paste0("Log2 Fold Change (Up in ", grp, ")")
  
  gg_bar <- ggplot(plot_data, aes(x = gene_name, y = log2FoldChange, fill = Type)) +
    geom_col(width = 0.7) + 
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("Significant (P < 0.05)" = "#CC3333", 
                                 "Not Significant" = "gray70")) +
    labs(title = plot_title, x = NULL, y = plot_y_label, fill = "Significance") +
    theme_bw() + 
    theme(
      axis.text.y = element_text(size = 11, color = "black", face = "bold.italic"),
      axis.text.x = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )
  
  pdf_file_name <- paste0("Barplot_Target_", grp, "_vs_WT_NoFilter.pdf")
  ggsave(pdf_file_name, plot = gg_bar, width = 8, height = 8)
  
  cat("  -> 已完成组别:", grp, "\n")
}

cat("\n=== 所有 Col Group 样本的靶基因分析完毕！ ===\n")

# ------------------------------------------------------------------------------
# 6. 计算 VIM2 和 VIM4 的联合表达水平 (仅保留的 OX 样本)
# ------------------------------------------------------------------------------
cat("\n=== 计算 VIM2+VIM4 合并表达量 (Col OX) ===\n")

norm_counts <- counts(dds_group, normalized = TRUE)
wt_samples <- rownames(colData_group)[colData_group$Group == "Control_Col_WT"]

vim_combined_wt <- norm_counts["AT1G66050", wt_samples] + norm_counts["AT1G66040", wt_samples]
wt_mean_combined <- mean(vim_combined_wt, na.rm = TRUE)

ox_samples <- rownames(colData_group)[colData_group$Group != "Control_Col_WT"]

combined_vim_results <- data.frame()

for (single_ox in ox_samples) {
  ox_combined_val <- norm_counts["AT1G66050", single_ox] + norm_counts["AT1G66040", single_ox]
  fold_change <- ox_combined_val / wt_mean_combined
  
  combined_vim_results <- rbind(combined_vim_results, data.frame(
    Sample_Name = single_ox,
    Group_Type = as.character(colData_group[single_ox, "Group"]),
    WT_Combined_Mean = wt_mean_combined,
    Sample_Combined_Value = ox_combined_val,
    Combined_Log2FC = log2(fold_change)
  ))
}

output_filename <- "Combined_VIM2_VIM4_Expression_and_Log2FC_All_OX_Grouped.txt"
write.table(combined_vim_results, file = output_filename, quote = FALSE, sep = "\t", row.names = FALSE)

cat("  -> 联合计算完成！结果已导出至:", output_filename, "\n")