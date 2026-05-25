
# ==============================================================================
# Script Name: plot_Figure5E.R
# Description: Generates Figure 5E - divergence_ATAC_data
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

##new fig5e ATAC divergence plot

rm(list=ls())



library(ggplot2)
library(data.table)

## Loading source code


source("/mnt/int/RIL/for_paper/codes_github/scripts/plotFits-overlay_interceptZero_formal_2021-09_using_genotypes_than_stress.R")

setwd("/mnt/int/RIL/for_paper/codes_github/data/")
input.data.dir<-"/mnt/int/RIL/for_paper/codes_github/data/"
output.data.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
###col
#PLEASE NAME YOUR FILES SUCH THAT YOU HAVE THE "name" for e.g Col-Mock followed by "category" for e.g global
global <- c("RIL50_ATAC_ABneutral_CG_estimates.Rdata",
            "RIL149_ATAC_ABneutral_CG_estimates.Rdata")


pedigree.list <- data.frame(global)


test_plotFITS(pedigree.names = pedigree.list,
              input.dir=input.data.dir,
              output.dir=output.data.dir,
              plot.type="both", #options for plot.type = "fit.only", "data.only", "both"
              out.name="Fig5E_divergence_ATAC_data",
              lsq.line="theory", #options for lsq.line = "theory", "pred" (don't change)
              alpha=0.3,
              facets.ncols=1, #because you have 4 categories: global, genes, TEs, promoters
              geom.point.size=1.6,
              geom.line.size=0.6,
              plot.height=4,
              plot.width=5)




