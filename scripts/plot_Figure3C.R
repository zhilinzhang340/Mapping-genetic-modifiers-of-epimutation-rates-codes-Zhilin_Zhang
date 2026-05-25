

# ==============================================================================
# Script Name: plot_Figure3C.R
# Description: Generates Figure 3C - gain rate loss rate lod curve for 68RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#Fig3C alpha beta lod curve add intergenic teM 2025-04 


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

data_file<- "newrun_2023-11_68_RIL_all_annotations_add_intergeneic_alpha_rate_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2025-04.csv"
my_cross <- read.cross(file = data_file, format = "csv", genotypes = c("0", "2"))
summary(my_cross)

my_cross_est <- est.map(my_cross,error.prob=1e-07)

my_cross <-replace.map(my_cross, my_cross_est)
my_cross

head(my_cross$pheno)


my_cross<- calc.genoprob(my_cross, step=1)

scan_global <- scanone(my_cross, pheno.col = "global", method = "hk")
scan_gene <- scanone(my_cross, pheno.col = "gene", method = "hk")
scan_TE <- scanone(my_cross, pheno.col = "TE", method = "hk")
scan_gbM <- scanone(my_cross, pheno.col = "gbM", method = "hk")
scan_UM <- scanone(my_cross, pheno.col = "UM", method = "hk")


scan_global$category <- "global"
scan_gene$category <- "gene"
scan_TE$category <- "TE"
scan_gbM$category <- "gbM"
scan_UM$category <- "UM"


scan_all <-rbind(scan_gene,scan_TE,scan_gbM,scan_UM)

perm_results_global <- scanone(my_cross, pheno.col = "global", n.perm = 1000)
su_perm_results_global <-summary(perm_results_global)
su_perm_results_global
cutoff_global <- su_perm_results_global[1,1]
cutoff_global

perm_results_gene <- scanone(my_cross, pheno.col = "gene", n.perm = 1000)
su_perm_results_gene <-summary(perm_results_gene)
su_perm_results_gene
cutoff_gene <- su_perm_results_gene[1,1]


perm_results_TE <- scanone(my_cross, pheno.col = "TE", n.perm = 1000)
su_perm_results_TE <-summary(perm_results_TE)
su_perm_results_TE
cutoff_TE <- su_perm_results_TE[1,1]

perm_results_gbM <- scanone(my_cross, pheno.col = "gbM", n.perm = 1000)
su_perm_results_gbM <-summary(perm_results_gbM)
su_perm_results_gbM
cutoff_gbM <- su_perm_results_gbM[1,1]

perm_results_UM <- scanone(my_cross, pheno.col = "UM", n.perm = 1000)
su_perm_results_UM <-summary(perm_results_UM)
su_perm_results_UM
cutoff_UM <- su_perm_results_UM[1,1]


library(ggplot2)
library(scales)




lodint(scan_global, chr=1, drop=2)
# chr      pos      lod category
# c1.loc163   1 163.0000 1.582873   global ->1_2346 
# 1_2623      1 184.9254 3.615377   global
# c1.loc190   1 190.0000 1.598616   global ->1_2661
scan_global[250:350,]

max(scan_global$lod)


# Step 1: 创建 QTL 对象
qtl_global <- makeqtl(cross = my_cross,
                      chr = 1,
                      pos = 184.9254,
                      what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "global"
fit_result_global <- fitqtl(cross = my_cross,
                            qtl = qtl_global,
                            pheno.col = "global",
                            method = "hk")

# Step 3: 输出结果
summary(fit_result_global)



scan_globalm <- scan_global[scan_global$lod==max(scan_global$lod),]
scan_globalm$annotate_label <- paste0("lod:",round(scan_globalm$lod,2),"\nInt:21.6-25.5Mb")
scan_globalm 



lodint(scan_gene, chr=1, drop=2)
# chr      pos      lod category
# c1.loc167   1 167.0000 2.474522     gene ->1_2417
# 1_2623      1 184.9254 4.585884     gene
# 1_2634      1 187.9566 2.301323     gene


# Step 1: 创建 QTL 对象
qtl_gene <- makeqtl(cross = my_cross,
                    chr = 1,
                    pos = 184.9254,
                    what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "gene"
fit_result_gene <- fitqtl(cross = my_cross,
                          qtl = qtl_gene,
                          pheno.col = "gene",
                          method = "hk")

# Step 3: 输出结果
summary(fit_result_gene)

scan_gene[250:350,]

post <- fread("/mnt/int/RIL/for_paper/codes_github/data/newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt")[,1:4]
start <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2417",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2634",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_gene$lod)

scan_genem <- scan_gene[scan_gene$lod==max(scan_gene$lod),]
scan_genem$annotate_label <- paste0("lod:",round(scan_genem$lod,2),"\nInt:",start,"-",end,"Mb")
scan_genem 


lodint(scan_gbM, chr=1, drop=2)
# chr      pos      lod category
# 1_2599      1 181.9401 4.619570      gbM
# 1_2623      1 184.9254 6.649192      gbM
# c1.loc187   1 187.0000 4.620353      gbM ->1_2634


# Step 1: 创建 QTL 对象
qtl_gbM <- makeqtl(cross = my_cross,
                   chr = 1,
                   pos = 184.9254,
                   what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "gbM"
fit_result_gbM <- fitqtl(cross = my_cross,
                         qtl = qtl_gbM,
                         pheno.col = "gbM",
                         method = "hk")

# Step 3: 输出结果
summary(fit_result_gbM)


scan_gbM[250:350,]
start_gbM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2599",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_gbM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2634",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_gbM$lod)


scan_gbMm <- scan_gbM[scan_gbM$lod==max(scan_gbM$lod),]
scan_gbMm$annotate_label <- paste0("lod:",round(scan_gbMm$lod,2),"\nInt:",start_gbM,"-",end_gbM,"Mb")
scan_gbMm 


scan_teM <- scanone(my_cross, pheno.col = "teM", method = "hk")
scan_teM$category <- "teM"
perm_results_teM <- scanone(my_cross, pheno.col = "teM", n.perm = 1000)
su_perm_results_teM <-summary(perm_results_teM)
su_perm_results_teM
cutoff_teM <- su_perm_results_teM[1,1]



scan_intergenic <- scanone(my_cross, pheno.col = "intergenic", method = "hk")
scan_intergenic$category <- "intergenic"
perm_results_intergenic <- scanone(my_cross, pheno.col = "intergenic", n.perm = 1000)
su_perm_results_intergenic <-summary(perm_results_intergenic)
su_perm_results_intergenic
cutoff_intergenic <- su_perm_results_intergenic[1,1]

#global +annotations

scan_all <-rbind(scan_global,scan_gene,scan_TE,scan_gbM,scan_UM,scan_teM,scan_intergenic)


scan_all$category <- factor(scan_all$category,c("global","gene","gbM","UM","TE","teM","intergenic"))


output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 10
plot.width <- 25
out.name <-"Fig3C_alpha_rate_global_add_annotations_add_intergeneic_teM_lod curves"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(scan_all, aes(x=pos, y=lod)) +  
  geom_line(aes(color=category))+
  scale_color_manual(values=c("black","#663300","blue","purple","pink","darkgreen","grey"))+
  labs(title="global & annotations alpha rate") +
  geom_hline(yintercept = cutoff_global, linetype = "dashed", color = "black", size = 0.2) +
  geom_text(data =scan_globalm, aes(label = annotate_label), x =scan_globalm$pos-15, y = scan_globalm$lod+0.3, hjust = 1, vjust =1, color = "black") +
  geom_hline(yintercept = cutoff_gene, linetype = "dashed", color = "#663300", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_TE, linetype = "dashed", color = "pink", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_gbM, linetype = "dashed", color ="blue", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_UM, linetype = "dashed", color = "purple", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_teM, linetype = "dashed", color = "darkgreen", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_intergenic, linetype = "dashed", color = "grey", size = 0.2) +  # 添加水平线
  
  geom_text(data =scan_genem, aes(label = annotate_label), x =scan_genem$pos-10, y = scan_genem$lod+0.3, hjust = 1, vjust =1, color ="#663300") +
  geom_text(data =scan_gbMm, aes(label = annotate_label), x =scan_gbMm$pos-10, y = scan_gbMm$lod+0.3, hjust = 1, vjust =1, color ="blue") +
  theme_classic() +
  facet_wrap(~ chr, nrow=1, scales = 'free_x') +
  scale_y_continuous(expand = c(0,0),limits = c(0,15)) +  # 使用 label_comma 避免科学计数法
  theme(text = element_text(size = 15),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )   # 添加次要网格线
print(gg)
dev.off()





#beta
rm(list=ls())
library(qtl)
library(qtlcharts)
library(data.table)
library(ggplot2)
library(stringr)
library(scales)  # 如果还没有加载 scales 包，需要先加载


post <- fread("/mnt/int/RIL/for_paper/codes_github/data/newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt")[,1:4]

#error.prob=1e-07

#run qtl

setwd("/mnt/int/RIL/for_paper/codes_github/data")

data_file<-"newrun_2023-11_68_RIL_all_annotations_add_intergeneic_beta_rate_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2025-04.csv"
my_cross <- read.cross(file = data_file, format = "csv", genotypes = c("0", "2"))
summary(my_cross)




my_cross_est <- est.map(my_cross,error.prob=1e-07)

my_cross <-replace.map(my_cross, my_cross_est)
my_cross

summary.map(my_cross)



my_cross<- calc.genoprob(my_cross, step=1)

scan_global <- scanone(my_cross, pheno.col = "global", method = "hk")
scan_gene <- scanone(my_cross, pheno.col = "gene", method = "hk")
scan_TE <- scanone(my_cross, pheno.col = "TE", method = "hk")
scan_gbM <- scanone(my_cross, pheno.col = "gbM", method = "hk")
scan_UM <- scanone(my_cross, pheno.col = "UM", method = "hk")


scan_global$category <- "global"
scan_gene$category <- "gene"
scan_TE$category <- "TE"
scan_gbM$category <- "gbM"
scan_UM$category <- "UM"


perm_results_global <- scanone(my_cross, pheno.col = "global", n.perm = 1000)
su_perm_results_global <-summary(perm_results_global)
su_perm_results_global
cutoff_global <- su_perm_results_global[1,1]
cutoff_global

perm_results_gene <- scanone(my_cross, pheno.col = "gene", n.perm = 1000)
su_perm_results_gene <-summary(perm_results_gene)
su_perm_results_gene
cutoff_gene <- su_perm_results_gene[1,1]


perm_results_TE <- scanone(my_cross, pheno.col = "TE", n.perm = 1000)
su_perm_results_TE <-summary(perm_results_TE)
su_perm_results_TE
cutoff_TE <- su_perm_results_TE[1,1]

perm_results_gbM <- scanone(my_cross, pheno.col = "gbM", n.perm = 1000)
su_perm_results_gbM <-summary(perm_results_gbM)
su_perm_results_gbM
cutoff_gbM <- su_perm_results_gbM[1,1]

perm_results_UM <- scanone(my_cross, pheno.col = "UM", n.perm = 1000)
su_perm_results_UM <-summary(perm_results_UM)
su_perm_results_UM
cutoff_UM <- su_perm_results_UM[1,1]



scan_all <-rbind(scan_gene,scan_TE,scan_gbM,scan_UM)

library(ggplot2)
library(scales)


lodint(scan_global, chr=1, drop=2)
# chr      pos      lod category
# c1.loc171   1 171.0000 4.038725   global ->1_2447
# 1_2623      1 184.9254 6.088450   global
# 1_2634      1 187.9566 4.003222   global
scan_global[250:350,]

post <- fread("/mnt/int/RIL/for_paper/codes_github/data/newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt")[,1:4]
start_global <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2447",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_global <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2634",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_global$lod)



# Step 1: 创建 QTL 对象
qtl_global <- makeqtl(cross = my_cross,
                      chr = 1,
                      pos = 184.9254,
                      what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "global"
fit_result_global <- fitqtl(cross = my_cross,
                            qtl = qtl_global,
                            pheno.col = "global",
                            method = "hk")

# Step 3: 输出结果
summary(fit_result_global)









scan_globalm <- scan_global[scan_global$lod==max(scan_global$lod),]
scan_globalm$annotate_label <- paste0("lod:",round(scan_globalm$lod,2),"\nInt:",start_global,"-",end_global,"Mb")
scan_globalm 



lodint(scan_gene, chr=1, drop=2)
# chr      pos       lod category
# 1_2601      1 183.4327  8.160926     gene
# 1_2623      1 184.9254 11.084320     gene
# c1.loc187   1 187.0000  8.777438     gene -> 1_2634
scan_gene[250:350,]
start_gene <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2601",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_gene <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2634",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_gene$lod)



# Step 1: 创建 QTL 对象
qtl_gene <- makeqtl(cross = my_cross,
                    chr = 1,
                    pos = 184.9254,
                    what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "gene"
fit_result_gene <- fitqtl(cross = my_cross,
                          qtl = qtl_gene,
                          pheno.col = "gene",
                          method = "hk")

# Step 3: 输出结果
summary(fit_result_gene)




scan_genem <- scan_gene[scan_gene$lod==max(scan_gene$lod),]
scan_genem$annotate_label <- paste0("lod:",round(scan_genem$lod,2),"\nInt:",start_gene,"-",end_gene,"Mb")
scan_genem 


lodint(scan_gbM, chr=1, drop=2)
# chr      pos      lod category
# c1.loc184   1 184.0000 12.39172      gbM -> 1_2601
# 1_2623      1 184.9254 14.93167      gbM
# c1.loc187   1 187.0000 12.21889      gbM ->1_2634
scan_gbM[250:350,]
start_gbM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2601",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_gbM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2634",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_gbM$lod)



# Step 1: 创建 QTL 对象
qtl_gbM <- makeqtl(cross = my_cross,
                   chr = 1,
                   pos = 184.9254,
                   what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "gbM"
fit_result_gbM <- fitqtl(cross = my_cross,
                         qtl = qtl_gbM,
                         pheno.col = "gbM",
                         method = "hk")

# Step 3: 输出结果
summary(fit_result_gbM)






scan_gbMm <- scan_gbM[scan_gbM$lod==max(scan_gbM$lod),]
scan_gbMm$annotate_label <- paste0("lod:",round(scan_gbMm$lod,2),"\nInt:",start_gbM,"-",end_gbM,"Mb")
scan_gbMm 

lodint(scan_TE, chr=1, drop=2)
# chr      pos       lod category
# c1.loc59    1  59.0000 0.8764704       TE. -> 1_581
# 1_2655      1 189.4493 2.9114659       TE
# c1.loc228   1 228.0000 0.7622921       TE. -> 1_2976


scan_TE[400:500,]
start_TE <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_581",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_TE <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2976",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_TE$lod)



scan_TEm <- scan_TE[scan_TE$lod==max(scan_TE$lod),]
scan_TEm$annotate_label <- paste0("lod:",round(scan_TEm$lod,2),"\nInt:",start_TE,"-",end_TE,"Mb")
scan_TEm 


# Step 1: 创建 QTL 对象
qtl_TE <- makeqtl(cross = my_cross,
                  chr = 1,
                  pos = 189.4493,
                  what = "prob")

# Step 2: 拟合模型，注意这里 pheno.col 要匹配新的列名 "TE"
fit_result_TE <- fitqtl(cross = my_cross,
                        qtl = qtl_TE,
                        pheno.col = "TE",
                        method = "hk")

# Step 3: 输出结果
summary(fit_result_TE)



scan_teM <- scanone(my_cross, pheno.col = "teM", method = "hk")
scan_teM$category <- "teM"
perm_results_teM <- scanone(my_cross, pheno.col = "teM", n.perm = 1000)
su_perm_results_teM <-summary(perm_results_teM)
su_perm_results_teM
cutoff_teM <- su_perm_results_teM[1,1]
cutoff_teM



scan_intergenic <- scanone(my_cross, pheno.col = "intergenic", method = "hk")
scan_intergenic$category <- "intergenic"
perm_results_intergenic <- scanone(my_cross, pheno.col = "intergenic", n.perm = 1000)
su_perm_results_intergenic <-summary(perm_results_intergenic)
su_perm_results_intergenic
cutoff_intergenic <- su_perm_results_intergenic[1,1]


#global +annotations

scan_all <-rbind(scan_global,scan_gene,scan_TE,scan_gbM,scan_UM,scan_teM,scan_intergenic)


scan_all$category <- factor(scan_all$category,c("global","gene","gbM","UM","TE","teM","intergenic"))



output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 10
plot.width <- 25
out.name <-"Fig3C_beta_rate_global_add_annotations_add_intergeneic_teM_lod curves"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(scan_all, aes(x=pos, y=lod)) +  
  geom_line(aes(color=category))+
  scale_color_manual(values=c("black","#663300","blue","purple","pink","darkgreen","grey"))+
  labs(title="global & annotations beta rate") +
  geom_hline(yintercept = cutoff_global, linetype = "dashed", color = "black", size = 0.2) +
  geom_text(data =scan_globalm, aes(label = annotate_label), x =scan_globalm$pos-15, y = scan_globalm$lod+0.3, hjust = 1, vjust =1, color = "black") +
  geom_hline(yintercept = cutoff_gene, linetype = "dashed", color ="#663300", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_TE, linetype = "dashed", color ="pink", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_gbM, linetype = "dashed", color = "blue", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_UM, linetype = "dashed", color ="purple", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_teM, linetype = "dashed", color = "darkgreen", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_intergenic, linetype = "dashed", color = "grey", size = 0.2) +  # 添加水平线
  
  geom_text(data =scan_genem, aes(label = annotate_label), x =scan_genem$pos-10, y = scan_genem$lod-0.3, hjust = 1, vjust =1, color ="#663300") +
  geom_text(data =scan_gbMm, aes(label = annotate_label), x =scan_gbMm$pos-10, y = scan_gbMm$lod-0.3, hjust = 1, vjust =1, color ="blue") +
  geom_text(data =scan_TEm, aes(label = annotate_label), x =scan_TEm$pos-30, y = scan_TEm$lod+0.7, hjust = 1, vjust =1, color ="pink") +
  theme_classic() +
  facet_wrap(~ chr, nrow=1, scales = 'free_x') +
  scale_y_continuous(expand = c(0,0),limits = c(0,15)) +  # 使用 label_comma 避免科学计数法
  theme(text = element_text(size = 15),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )   # 添加次要网格线
print(gg)
dev.off()



