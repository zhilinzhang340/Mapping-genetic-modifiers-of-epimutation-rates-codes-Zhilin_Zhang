
# ==============================================================================
# Script Name: plot_Figure4H.R
# Description: Generates Figure 4H - barplot Gene expression comparison between peak genotypes
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

rm(list=ls())



library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
combined_df2 <- fread("/mnt/int/RIL/for_paper/codes_github/data/preTableSNew_expression_level_across_67RILs.txt")



combined_df20 <- combined_df2[combined_df2$peak_genotype=="peakCvi",]
combined_df22 <- combined_df2[combined_df2$peak_genotype=="peakLer",]


usesler <- combined_df22[,-c(1,21)]

usescvi <- combined_df20[,-c(1,21)]

usesler$genotype <- "peakLer"
usescvi$genotype <- "peakCvi"

usesall <- rbind(usesler,usescvi)
head(usesall )






# Step 1: 数据转换为 long format
#------------------------
usesall_long <- usesall %>%
  pivot_longer(cols = -genotype, names_to = "Gene", values_to = "Expression")

#------------------------
# Step 2: 计算每个 Gene 的 t 检验 p 值
#------------------------
p_values <- compare_means(Expression ~ genotype,
                          group.by = "Gene",
                          data = usesall_long,
                          method = "t.test")

#------------------------
# Step 3: 标记显著性 TRUE/FALSE
#------------------------
p_values <- p_values %>%
  mutate(is_significant = p < 0.05)

#------------------------
# Step 4: 生成每个基因的显著性标注位置
#------------------------
y_max <- usesall_long %>%
  group_by(Gene) %>%
  summarise(y_pos = max(Expression) + 0.5)

# 显著性标注数据
p_sig <- p_values %>%
  filter(is_significant) %>%
  left_join(y_max, by = "Gene")

#------------------------
# Step 5: 合并显著性标签回原始表达数据
#------------------------
usesall_long <- usesall_long %>%
  left_join(p_values %>% dplyr::select(Gene, is_significant), by = "Gene")

col_names_ordered <-row_names_ordered <- c("Vim4","MET1","Vim1","CMT3","Vim3","Vim5","Vim6","CMT1","CMT2","DRM1","DRM2",
                                           "DML2","DML3","DME","ROS1","ROS3","ROS4","ROS5","REF6")


usesall_long$Gene <- factor(usesall_long$Gene,levels=col_names_ordered )
usesall_long$genotype <- factor(usesall_long$genotype,levels=c("peakLer","peakCvi"))
# 设置 baseline
baseline <- 7

# Step A: 将表达量减去 baseline
usesall_long <- usesall_long %>%
  mutate(Expression_adj = Expression - baseline)

# Step B: 调整显著性标注位置（同样以 baseline 为基准）
y_max <- usesall_long %>%
  group_by(Gene) %>%
  summarise(y_pos = max(Expression_adj) + 0.2)

p_sig <- p_values %>%
  filter(is_significant) %>%
  left_join(y_max, by = "Gene")


output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 8
out.name<-"Fig4H barplot Gene expression comparison between peak genotypes"
pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

# Step C: 绘图
gg <- ggplot(usesall_long, aes(x = Gene, y = Expression_adj, fill = genotype, alpha = is_significant)) +
  geom_bar(stat = "summary", fun = mean, position = position_dodge(0.8), width = 0.7) +
  
  geom_errorbar(stat = "summary", fun.data = mean_se,
                position = position_dodge(0.8),
                width = 0.2, linewidth = 0.3) +
  
  geom_text(data = p_sig,
            aes(x = Gene, y = y_pos, label = p.signif),
            inherit.aes = FALSE,
            vjust = 5, size = 5) +
  
  # ✅ y轴从0开始，但显示原始表达值（7起）
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 4),  # 可根据最大值微调
    breaks = seq(0, 4, by = 0.5),
    labels = seq(baseline, baseline + 4, by = 0.5)
  ) +
  
  scale_fill_manual(values = c("peakLer" = "#1D2088", "peakCvi" = "#00A29A")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
  
  theme_classic() +
  labs(
    title = "Gene expression comparison between genotypes",
    x = "",
    y = "Expression level (mean ± SE)"
  ) +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


print(gg)
dev.off()
