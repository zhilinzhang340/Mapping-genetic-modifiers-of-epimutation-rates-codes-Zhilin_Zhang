

# ==============================================================================
# Script Name: plot_Figure5G.R
# Description: Generates Figure 5G - phenotype with peak markers 68 RIL 312 traits
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


rm(list=ls())
gc()

library(ggplot2)
p_value_df <- fread("/mnt/int/RIL/for_paper/codes_github/data/phenotype_QTL_68_RIL_312_traits.txt")

# 使用p.adjust函数来计算FDR
adjusted_p_values <- p.adjust(p_value_df$P_Value, method = "BH")

# 将调整后的p值添加到数据框中
p_value_df$adjusted_p_values <- adjusted_p_values
p_value_df[p_value_df$adjusted_p_values<0.05,]


# 将阈值转换为相应的-log10(FDR)值，以便在图表中使用
fdr_log_threshold <- -log10(min(adjusted_p_values[adjusted_p_values < 0.05], na.rm = TRUE))

# 打印FDR阈值
print(fdr_log_threshold)


head(p_value_df)
p_value_df$logadj <- -log10(p_value_df$adjusted_p_values)


output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 10
plot.width <- 20
out.name<-"Fig5G phenotype with peak markers 68 RIL 312 traits adjusted_p_values"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(data = p_value_df, mapping = aes(x = trait, y = logadj)) +
  geom_col(width = .7, position = 'dodge') +
  theme_classic() +
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  labs(y = "-log10(adjusted_p_values)", x = "", title = "312 traits with peak marker") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, max(p_value_df$trait, na.rm = TRUE), by = 5),
                     expand = c(0, 0), 
                     limits = c(0,max(p_value_df$trait)+5))+
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 0.8) + # 添加红色虚线
  annotate("text", x = 20, y = -log10(0.05) + 0.05, 
           label = "adjusted_p_values = 0.05", color = "red", size = 4)  # 添加文本注释


print(gg)
dev.off()


