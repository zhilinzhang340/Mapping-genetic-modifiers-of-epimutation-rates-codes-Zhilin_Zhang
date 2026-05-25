

#Fig1C
# ==============================================================================
# Script Name: plot_Figure1C.R
# Description: Generates Figure 1C - divergence plot for MA-accessions.
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================
##formal
##CG Col
rm(list=ls())

library(ggplot2)
library(data.table)

## Loading source code

#source("~/documents/MA-lines/Rscript/plotFits-overlay_interceptZero_formal_2021-09.R")
source("/mnt/int/RIL/for_paper/codes_github/scripts/plotFits-overlay_interceptZero_formal_2021-09_using_genotypes_than_stress_adjust_2023-11.R")

setwd("/mnt/int/RIL/for_paper/codes_github/data/")
input.data.dir<-"/mnt/int/RIL/for_paper/codes_github/data/"
output.data.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"

global <- c("Col-0_Site_global_ABneutral_CG_estimates.Rdata",
            "Kn-0_Site_global_ABneutral_CG_estimates.Rdata",
            "Mt-0_Site_global_ABneutral_CG_estimates.Rdata",
            "Tsu-0_Site_global_ABneutral_CG_estimates.Rdata")

pedigree.list <- data.frame(global)


test_plotFITS(pedigree.names = pedigree.list,
              input.dir=input.data.dir,
              output.dir=output.data.dir,
              plot.type="both", #options for plot.type = "fit.only", "data.only", "both"
              out.name="Fig1C divergence_4_accessions_site_only_global",
              lsq.line="theory", #options for lsq.line = "theory", "pred" (don't change)
              alpha=0.3,
              facets.ncols=1, #because you have 4 categories: global, genes, TEs, promoters
              geom.point.size=1,
              geom.line.size=0.6,
              plot.height=4,
              plot.width=6)
