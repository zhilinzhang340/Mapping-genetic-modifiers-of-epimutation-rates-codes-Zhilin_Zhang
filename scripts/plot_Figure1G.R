

# ==============================================================================
# Script Name: plot_Figure1G.R
# Description: Generates Figure 1G - old and new genetic map.
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


#old genetic map

rm(list=ls())


library(qtl)
library(qtlcharts)
setwd("/mnt/int/RIL/for_paper/codes_github/data/")

data_file<- "/mnt/int/RIL/for_paper/codes_github/data/known_makers_all_Beta_rate_add_genotype_matrix_without_ID_only_68RIL_lines_all_annotations_2023-11.csv"
my_cross <- read.cross(file = data_file, format = "csv", genotypes = c("A", "B"))
summary(my_cross)


my_cross_est <- est.map(my_cross,error.prob=5e-04)

my_cross <-replace.map(my_cross, my_cross_est)


output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 8
plot.width <- 8
out.name <-"Fig1G old genetic map 144 known markers"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- plot.map(my_cross)

print(gg)
dev.off()




# new genetic map
rm(list=ls())
library(qtl)
library(qtlcharts)
library(data.table)
library(ggplot2)
library(stringr)
library(scales)  # 如果还没有加载 scales 包，需要先加载

#error.prob=1e-07

#run qtl

setwd("/mnt/int/RIL/for_paper/codes_github/data/")

data_file<- "newrun_2023-11_68_RIL_all_annotations_alpha_rate_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2023-11.csv"
my_cross <- read.cross(file = data_file, format = "csv", genotypes = c("0", "2"))
summary(my_cross)

my_cross_est <- est.map(my_cross,error.prob=1e-07)

my_cross <-replace.map(my_cross, my_cross_est)
my_cross


out_map <- summary.map(my_cross)

output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 8
plot.width <- 8
out.name <-"Fig1G new genetic map 732 SNP markers"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- plot.map(my_cross)

print(gg)
dev.off()



