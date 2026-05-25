
# ==============================================================================
# Script Name: plot_Figure5F.R
# Description: Generates Figure 5F - Mean Variance_with_SE all 312 traits
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================
rm(list=ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)
outs <-read.csv("/mnt/int/RIL/for_paper/codes_github/data/pheno_peak_status_68RIL_312phenos_prepare_lm_1_2623_0as-1_2as1.csv")
outs

# 替换"-"为NA（仅处理第三列及之后的所有列）
outs[, 3:ncol(outs)] <- lapply(outs[, 3:ncol(outs)], function(x) replace(x, x == "-", NA))
outs[, 3:ncol(outs)] <- lapply(outs[, 3:ncol(outs)], as.numeric)  # 确保转换为数值
use <- outs

use$type <- ifelse(use$peak_status == "-1", "Cvi",
                   ifelse(use$peak_status == "1", "Ler", NA))

# 初始化存储结果
results <- data.frame(Trait = character(), Group = character(), Variance = numeric(), stringsAsFactors = FALSE)

# 遍历所有 traits（第3列到最后）
for (i in 3:(ncol(use)-1)) {

  trait_name <- colnames(use)[i]
  cvi_vals <- use[use$type == "Cvi", i]
  ler_vals <- use[use$type == "Ler", i]
  
  cvi_var <- var(cvi_vals, na.rm = TRUE)
  ler_var <- var(ler_vals, na.rm = TRUE)

  
  results <- rbind(results,
                   data.frame(Trait = trait_name, Group = "Cvi", Variance = cvi_var),
                   data.frame(Trait = trait_name, Group = "Ler", Variance = ler_var))
}




results2 <- results[!is.na(results$Variance), ]




# STEP 1: 数据清理并转换为宽格式，准备 paired t-test
results2_clean <- results2 %>%
  group_by(Trait) %>%
  filter(n() == 2, all(!is.na(Variance))) %>%
  ungroup()

results_wide <- pivot_wider(results2_clean, names_from = Group, values_from = Variance)

# STEP 2: 执行 paired t-test
t_test_result <- t.test(results_wide$Cvi, results_wide$Ler, paired = TRUE)
pval <- t_test_result$p.value
p_text <- ifelse(pval < 0.001, "p < 0.001",
                 ifelse(pval < 0.01, "p < 0.01",
                        ifelse(pval < 0.05, "p < 0.05", paste0("p = ", round(pval, 3)))))

# STEP 3: 计算绘图所需的 mean 和 SE
summary_stats <- results2_clean %>%
  group_by(Group) %>%
  summarise(
    mean_variance = mean(Variance, na.rm = TRUE),
    se = sd(Variance, na.rm = TRUE) / sqrt(n())
  )





# STEP 4: 设置绘图参数
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
out.name <- "Fig5F Mean Variance_with_SE all 312 traits"
plot.height <- 5
plot.width <- 5

pdf(paste0(output.dir, out.name, ".pdf"), width = plot.width, height = plot.height)

# STEP 5: 绘图 + 添加显著性文字
gg <- ggplot(summary_stats, aes(x = Group, y = mean_variance, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_errorbar(aes(ymin = mean_variance - se, ymax = mean_variance + se),
                width = 0.2, position = position_dodge(0.6)) +
  labs(title = "Mean variance all 312 traits",
       x = "Genotype group", y = "Mean variance ± SE") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_classic() +
  geom_text(aes(label = format(round(mean_variance, -4), scientific = TRUE),
                y = mean_variance + se + 0.05 * max(mean_variance)),
            position = position_dodge(0.6),
            size = 5) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,11e+06)) +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),legend.position = "none") +
  # 添加显著性文字
  annotate("text", x = 1.5, y = max(summary_stats$mean_variance + summary_stats$se) * 1.1,
           label = p_text, size = 5)

print(gg)
dev.off()



