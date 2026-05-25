

# ==============================================================================
# Script Name: plot_Figure5B.R
# Description: Generates Figure 5B - cvi_accessions_var_methy_lvl_within_without_deletion
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================
rm(list=ls())
gc()

use <- fread("/mnt/int/RIL/for_paper/codes_github/data/TableSNew_Cvi_accssions_global_mCG_level_deletion.txt")


use$Deletion_status <- factor(use$Deletion_status,levels = c("Without_deletion","With_deletion"))
head(use)



use$global_mCG_level <- 100*use$global_mCG_level

summary_df <- use %>%
  group_by(Deletion_status) %>%
  summarise(
    n = n(),
    var_level = var(global_mCG_level, na.rm = TRUE),
    se_var = var_level * sqrt(2 / (n - 1)),
    se_lower = var_level - se_var,
    se_upper = var_level + se_var
  )

summary_df 

summary_df$Deletion_status <- factor(summary_df$Deletion_status,levels = c("Without_deletion","With_deletion")) 


var.test(global_mCG_level ~ Deletion_status, data = use)
ftest <- var.test(global_mCG_level ~ Deletion_status, data = use)
p_value <- ftest$p.value
p_value 

output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name <- "Fig5B_cvi_accessions_var_methy_lvl_within_without_deletion_with_SE"

pdf(paste0(output.dir, out.name, ".pdf"), width = plot.width, height = plot.height)

gg <- ggplot(summary_df, aes(x = Deletion_status, y = var_level, fill = Deletion_status)) +
  geom_bar(stat = "identity", width = 0.6,alpha=0.8) +
  geom_errorbar(aes(ymin = se_lower, ymax = se_upper), width = 0.2, size = 0.5) +
  geom_text(aes(y=se_upper+0.1,label = round(var_level, 2)), vjust = -0.5, size = 5) +
  labs(title = "83 Cvi accessions group by deletion",
       y = "Variance of methylation level across accessions", x = "Deletion between VIM2 and VIM4") +
  theme_classic() +
  scale_fill_manual(values = c("Without_deletion" = "#4DBBD5", "With_deletion" = "#E64B35")) +
  annotate("text", x = 1.5, y = 6.3, 
           label = paste0("F-test p = ", signif(p_value, 3)),
           size = 5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6.5)) +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

print(gg)
dev.off()

