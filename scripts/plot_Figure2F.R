

# ==============================================================================
# Script Name: plot_Figure2F.R
# Description: Generates Figure 2F - global_per5kb genetic Constitution Percentage
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================




### make gene TE annotation for trans qtl global window chr1 chr2

rm(list=ls())


library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(IRanges)
library(Biostrings)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales)

setwd("/mnt/int/RIL/for_paper/codes_github/data")

gbM<-fread("gbM_gene_anotation_extract_Arabidopsis.bed")
teM<-fread("teM_gene_anotation_extract_Arabidopsis.bed")
UM<-fread("UM_gene_anotation_extract_Arabidopsis.bed")
gene<-fread("gene_anotation_extract_Arabidopsis.bed")
other_gene<-fread("other_gene_execpt_gbM_teM_UM_TAIR10_Arabidopsis_whole_gene.bed")
TE<-fread("TAIR10_TE.bed")
intergenic <- import.gff3("intergenic.gff3",colnames=c("type", "ID","source"))
intergenic$type <- "intg"


colnames(gbM)<-c("seqnames","start","end","id","class","strand")
colnames(teM)<-c("seqnames","start","end","id","class","strand")
colnames(UM)<-c("seqnames","start","end","id","class","strand")
colnames(other_gene)<-c("seqnames","start","end","id","class","strand")
colnames(TE)<-c("seqnames","start","end","id","class","strand")


gr_gbM<-GRanges(gbM)
gr_teM<-GRanges(teM)
gr_UM<-GRanges(UM)
gr_other_gene<-GRanges(other_gene)
gr_TE<-GRanges(TE)



class<-c("gbM","teM","UM","other_gene","TE","intergenic","others")

df.all<-NULL


CS <- fread("68RIL_global_per5kb_qtl_peak_marker_chr1_with_trans_pos.txt")
CS <- CS[,13:15]

colnames(CS)<-c("seqnames","start","end")
gr<-GRanges(CS)

inter_gbM<-GenomicRanges::intersect(gr,gr_gbM,ignore.strand=T)
inter_teM<-GenomicRanges::intersect(gr,gr_teM,ignore.strand=T)
inter_UM<-GenomicRanges::intersect(gr,gr_UM,ignore.strand=T)
inter_other_gene<-GenomicRanges::intersect(gr,gr_other_gene,ignore.strand=T)
inter_TE<-GenomicRanges::intersect(gr,gr_TE,ignore.strand=T)
inter_intergenic<-GenomicRanges::intersect(gr,intergenic,ignore.strand=T)


percentage<-c(100*sum(width(inter_gbM))/sum(width(gr)),
              100*sum(width(inter_teM))/sum(width(gr)),
              100*sum(width(inter_UM))/sum(width(gr)),
              100*sum(width(inter_other_gene))/sum(width(gr)),

              100*sum(width(inter_TE))/sum(width(gr)), 

              100*sum(width(inter_intergenic))/sum(width(gr)),
              100*(1-sum(width(inter_gbM))/sum(width(gr))
                   -sum(width(inter_teM))/sum(width(gr))
                   -sum(width(inter_UM))/sum(width(gr))
                   -sum(width(inter_other_gene))/sum(width(gr))
                   -sum(width(inter_TE))/sum(width(gr))
                   -sum(width(inter_intergenic))/sum(width(gr))

              ))


df_5kbc1<-data.frame(class,"chr1",percentage,"5kb")





CS <- fread("68RIL_global_per5kb_qtl_peak_marker_chr2_with_trans_pos.txt")
CS <- CS[,13:15]


colnames(CS)<-c("seqnames","start","end")
gr<-GRanges(CS)

inter_gbM<-GenomicRanges::intersect(gr,gr_gbM,ignore.strand=T)
inter_teM<-GenomicRanges::intersect(gr,gr_teM,ignore.strand=T)
inter_UM<-GenomicRanges::intersect(gr,gr_UM,ignore.strand=T)
inter_other_gene<-GenomicRanges::intersect(gr,gr_other_gene,ignore.strand=T)
inter_TE<-GenomicRanges::intersect(gr,gr_TE,ignore.strand=T)
inter_intergenic<-GenomicRanges::intersect(gr,intergenic,ignore.strand=T)

percentage<-c(100*sum(width(inter_gbM))/sum(width(gr)),
              100*sum(width(inter_teM))/sum(width(gr)),
              100*sum(width(inter_UM))/sum(width(gr)),
              100*sum(width(inter_other_gene))/sum(width(gr)),

              100*sum(width(inter_TE))/sum(width(gr)), 

              100*sum(width(inter_intergenic))/sum(width(gr)),

              100*(1-sum(width(inter_gbM))/sum(width(gr))
                   -sum(width(inter_teM))/sum(width(gr))
                   -sum(width(inter_UM))/sum(width(gr))
                   -sum(width(inter_other_gene))/sum(width(gr))
                   -sum(width(inter_TE))/sum(width(gr))
                   -sum(width(inter_intergenic))/sum(width(gr))

              ))


df_5kbc2<-data.frame(class,"chr2",percentage,"5kb")




#ref
total <- fread("TAIR10_whole_chr_position.txt")


sumall <- sum(total$V2)

percentage<-c(100*sum(width(gr_gbM))/sumall,
              100*sum(width(gr_teM))/sumall,
              100*sum(width(gr_UM))/sumall,
              100*sum(width(gr_other_gene))/sumall,

              100*sum(width(gr_TE))/sumall, 

              100*sum(width(intergenic))/sumall, 

              100*(1-sum(width(gr_gbM))/sumall
                   -sum(width(gr_teM))/sumall
                   -sum(width(gr_UM))/sumall
                   -sum(width(gr_other_gene))/sumall
                   -sum(width(gr_TE))/sumall

                   -sum(width(intergenic))/sumall

              ))

df_ref<-data.frame(class,"all",percentage,"ref")

names(df_5kbc1) <- c("class","chr","percentage","window")
names(df_5kbc2) <- c("class","chr","percentage","window")
names(df_ref) <- c("class","chr","percentage","window")

df.all5kb<-rbind(df_5kbc1,df_5kbc2,df_ref)   
df.all5kb $class<-factor(df.all5kb $class,levels = c("gbM","teM","UM","other_gene","TE","intergenic","others"))



output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 6
plot.width <- 8
out.name<-"Fig2F global 5kb window chr1 chr2 trans genetic Constitution Percentage"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)



gg <- ggplot(df.all5kb, aes(x = chr, y = percentage, fill = class)) +
  
  geom_bar(stat = "identity", colour = "black") +
  
  labs( y ="Percentage",x="",title="global 5kb window chr1 chr2 trans") +
  
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  
  theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5,size=15))+
  geom_text(aes(label = round(percentage, 1)), 
            position = position_stack(vjust = 0.5), # 居中对齐
            size = 4) 


print(gg)
dev.off()

