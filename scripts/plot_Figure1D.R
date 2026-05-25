

#Fig1D
# ==============================================================================
# Script Name: plot_Figure1D.R
# Description: Generates Figure 1D - Global gain rate and loss rate of MA-accessions
# Author: Zhilin Zhang
# Date: 2026-01
# ==============================================================================
#rate

rm(list=ls())


library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales)


setwd("/mnt/int/RIL/for_paper/codes_github/data/")


#global
##Col
##alpha
global <- fread("/mnt/int/RIL/for_paper/codes_github/data/gain_loss-rate_MA-accessions_used_2026-01.txt")
globala <- global[order(global$Alpha),]

globala$Stress <- factor(globala$Stress ,globala$Stress)


output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name<-"Fig1D alpha rate 4 accessions CG global"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

gg<-ggplot(globala,aes(x=Stress,weight=Alpha/1e-04))+
  geom_bar( width = .7,position = 'dodge') +
  geom_errorbar(aes(ymin = Alpha/1e-04 - Alpha_SE/1e-04, ymax = Alpha/1e-04 + Alpha_SE/1e-04), width = 0.4, size = 0.3, position = position_dodge(0.7)) +
  labs( y =expression(paste('Gain Rate, ',alpha,' (x',10^-4,')')),x="",title="4 accessions CG global") +
  theme_classic()+
  geom_text(aes(y=(Alpha/1e-04+Alpha_SE/1e-04)+0.2,label=round(Alpha/1e-04,2)),size=3,position=position_dodge(0.7),vjust=-1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(globala$Alpha/1e-04 +globala$Alpha_SE/1e-04)+0.5))+
  theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5),)

print(gg)
dev.off()


globalb <- global[order(global$Beta),]

globalb$Stress <- factor(globalb$Stress ,globalb$Stress)

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name<-"Fig1D beta rate 4 accessions CG global"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

gg<-ggplot(globalb,aes(x=Stress,weight=Beta/1e-04))+
  geom_bar( width = .7,position = 'dodge') +
  geom_errorbar(aes(ymin = Beta/1e-04 - Beta_SE/1e-04, ymax = Beta/1e-04 + Beta_SE/1e-04), width = 0.4, size = 0.3, position = position_dodge(0.7)) +
  labs( y =expression(paste('Loss Rate, ',beta,' (x',10^-4,')')),x="",title="4 accessions CG global") +
  theme_classic()+
  geom_text(aes(y=(Beta/1e-04+Beta_SE/1e-04)+0.2,label=round(Beta/1e-04,2)),size=3,position=position_dodge(0.7),vjust=-1)+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(globalb$Beta/1e-04 +globalb$Beta_SE/1e-04)+1.5))+
  theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5))

print(gg)
dev.off()
