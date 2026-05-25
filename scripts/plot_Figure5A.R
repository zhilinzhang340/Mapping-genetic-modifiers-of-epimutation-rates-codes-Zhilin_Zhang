

# ==============================================================================
# Script Name: plot_Figure5A.R
# Description: Generates Figure 5A - boxplot cvi accessions mean methy lvl within without deletion
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

rm(list=ls())
gc()

use <- fread("/mnt/int/RIL/for_paper/codes_github/data/TableSNew_Cvi_accssions_global_mCG_level_deletion.txt")


use$Deletion_status <- factor(use$Deletion_status,levels = c("Without_deletion","With_deletion"))
head(use)


output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name <-"Fig5A boxplot cvi accessions mean methy lvl within without deletion"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

gg <-ggplot(use, aes(x = Deletion_status, y = 100 * global_mCG_level,color = Deletion_status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA) +  # 只画边框
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +  # 指定颜色
  #scale_fill_manual(values = c("Without_deletion" = "#4DBBD5", "With_deletion" = "#E64B35")) +
  scale_color_manual(values = c("Without_deletion" = "#4DBBD5", "With_deletion" = "#E64B35")) +
  labs(title = "83 Cvi accessions group by deletion",
       y = "Mean methylation level across accessions (%)", x = "Deletion between VIM2 and VIM4") +
  theme_classic() +
  scale_y_continuous(limits = c(13, 25)) +
  stat_compare_means(method = "t.test", label.y = 24) +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")


print(gg)
dev.off()

