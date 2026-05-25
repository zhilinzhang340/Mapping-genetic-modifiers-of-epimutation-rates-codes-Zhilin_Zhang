

# ==============================================================================
# Script Name: plot_Figure1A.R
# Description: Generates Figure 1A - Global mCG methylation level distribution 
#              of 850 accessions with parental lines marked.
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

# --- 1. Setup and Libraries ---
rm(list=ls()) # Clear environment

# Load necessary libraries
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
ma <- fread("/mnt/int/RIL/for_paper/codes_github/data/new_observed_methy_level_using_proportion_MA-accession_2023-11.txt")
ma$name <- paste0(str_split_fixed(ma$Genotypes,"0",2)[,1],"-0")


global <- fread("/mnt/int/RIL/for_paper/codes_github/data/RIL_371_all_sample_CG_new_methy_level_without_C_M_global_2026-01.txt")
globals <- global[global$line=="Cvi"|global$line=="Ler"]

# 1. 按 line 分组计算 Observed_methy_level 的平均值
rils <- globals[, .(mean_methy_lvl = mean(Observed_methy_level)), by = line]

# 2. 动态生成 name 列，例如将 Cvi 转换为 Cvi-0
rils[, name := paste0(line, '-0')]

# 3. 修改 line 列的格式，例如将 Cvi 转换为 Cvi0
rils[, line := paste0(line, '0')]

# 4. 调整列的顺序，使其与 rils 的列顺序完全一致
setcolorder(rils, c('line', 'mean_methy_lvl', 'name'))
rils




##850 samples methy level

data<-fread("/mnt/int/RIL/for_paper/codes_github/data/850_samples_CG_methy_level_used_2026-01.txt")



output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 4
plot.width <- 5 
out.name<-"Fig1A 850 samples CG methy level mark 6 accssions"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(data = data, aes(x = 100 * adjust)) +
  geom_histogram(bins = 30, alpha = 0.5) +
  geom_vline(data = rils, aes(xintercept = 100 * mean_methy_lvl), color = "red") +
  geom_vline(data = ma, aes(xintercept = 100 * MA_obs_lvl), color = "blue") +
  geom_text_repel(data = rils, aes(x = 100 *mean_methy_lvl,y =75, label = name,color = "red")) +
  geom_text_repel(data =ma, aes(x = 100 * MA_obs_lvl,y =55, label = name, color = "blue")) +
  labs(x = "mCG Methylation level global (%)", y = "count", title = "850 samples CG methy level") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),legend.position="none")+
  scale_color_manual(values = c("red" = "red",  "blue" = "blue")) 

print(gg)
dev.off()



