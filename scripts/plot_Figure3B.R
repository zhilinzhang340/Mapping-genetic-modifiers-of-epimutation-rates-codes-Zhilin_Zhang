

# ==============================================================================
# Script Name: plot_Figure3B.R
# Description: Generates Figure 3B - divergence plot RILs
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


##new RIL fig 3B divergence plot 2025-06 
##CG Col
rm(list=ls())

library(ggplot2)
library(data.table)

## Loading source code

source("/mnt/int/RIL/for_paper/codes_github/scripts/plotfits_divergence_adjust_fig3b_two_color_2025-05.R")

setwd("/mnt/int/RIL/for_paper/codes_github/data/")
input.data.dir<-"/mnt/int/RIL/for_paper/codes_github/data/"
output.data.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
###col
#PLEASE NAME YOUR FILES SUCH THAT YOU HAVE THE "name" for e.g Col-Mock followed by "category" for e.g global
use <- c( 58, 150, 182, 149,50,15,155,162)

uses <- paste0("RIL",use,"_Site_global_ABneutral_CG_estimates.Rdata")


pedigree.list <- data.frame(uses)
pedigree.list

test_plotFITS(pedigree.names = pedigree.list,
              input.dir=input.data.dir,
              output.dir=output.data.dir,
              plot.type="both", #options for plot.type = "fit.only", "data.only", "both"
              out.name="Fig3B divergence_4low_4high_site_only_global",
              lsq.line="theory", #options for lsq.line = "theory", "pred" (don't change)
              alpha=0.3,
              facets.ncols=1, #because you have 4 categories: global, genes, TEs, promoters
              geom.point.size=1.6,
              geom.line.size=0.6,
              plot.height=4,
              plot.width=5)

