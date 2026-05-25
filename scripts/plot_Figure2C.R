
#Fig2C
# ==============================================================================
# Script Name: plot_Figure2C.R
# Description: Generates Figure 2C - QTL LOD curve methylation level RILs
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#methy lvl all

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

data_file<-"newrun_2024-11_68_RIL_global_gbM_gene_UM_TE_intergenic_teM_new_methy_lvl_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2024-11.csv"
my_cross <- read.cross(file = data_file, format = "csv", genotypes = c("0", "2"))
summary(my_cross)


my_cross_est <- est.map(my_cross,error.prob=1e-07)

my_cross <-replace.map(my_cross, my_cross_est)
my_cross


my_cross<- calc.genoprob(my_cross, step=1)

scan_global <- scanone(my_cross, pheno.col = "global_methy_lvl", method = "hk")
scan_gbM <- scanone(my_cross, pheno.col = "gbM_methy_lvl", method = "hk")

scan_global$category <- "global_methy_lvl"
scan_gbM$category <- "gbM_methy_lvl"


perm_results_global <- scanone(my_cross, pheno.col = "global_methy_lvl", n.perm = 1000)
su_perm_results_global <-summary(perm_results_global)
su_perm_results_global
cutoff_global <- su_perm_results_global[1,1]


perm_results_gbM <- scanone(my_cross, pheno.col = "gbM_methy_lvl", n.perm = 1000)
su_perm_results_gbM <-summary(perm_results_gbM)
su_perm_results_gbM
cutoff_gbM <- su_perm_results_gbM[1,1]

scan_gene <- scanone(my_cross, pheno.col = "gene_methy_lvl", method = "hk")
scan_UM <- scanone(my_cross, pheno.col = "UM_methy_lvl", method = "hk")



scan_gene$category <- "gene_methy_lvl"
scan_UM$category <- "UM_methy_lvl"



perm_results_gene <- scanone(my_cross, pheno.col = "gene_methy_lvl", n.perm = 1000)
su_perm_results_gene <-summary(perm_results_gene)
su_perm_results_gene
cutoff_gene <- su_perm_results_gene[1,1]


perm_results_UM <- scanone(my_cross, pheno.col = "UM_methy_lvl", n.perm = 1000)
su_perm_results_UM <-summary(perm_results_UM)
su_perm_results_UM
cutoff_UM <- su_perm_results_UM[1,1]

scan_TE <- scanone(my_cross, pheno.col = "TE_methy_lvl", method = "hk")




scan_TE$category <- "TE_methy_lvl"




perm_results_TE <- scanone(my_cross, pheno.col = "TE_methy_lvl", n.perm = 1000)
su_perm_results_TE <-summary(perm_results_TE)
su_perm_results_TE
cutoff_TE <- su_perm_results_TE[1,1]



scan_intergenic <- scanone(my_cross, pheno.col = "intergenic_methy_lvl", method = "hk")



scan_intergenic$category <- "intergenic_methy_lvl"



perm_results_intergenic <- scanone(my_cross, pheno.col = "intergenic_methy_lvl", n.perm = 1000)
su_perm_results_intergenic <-summary(perm_results_intergenic)
su_perm_results_intergenic
cutoff_intergenic <- su_perm_results_intergenic[1,1]



lodint(scan_global, chr=1, drop=2)
# chr      pos      lod         category
# c1.loc179   1 179.0000 6.468837 global_methy_lvl -> 1_2543 select the cloest one before 
# 1_2661      1 190.9419 8.616885 global_methy_lvl
# 1_2690      1 196.9125 6.611722 global_methy_lvl 


#check major QTL was identified on chromosome 1, explaining  x% of total CG methylation level variance
library(qtl)

# Step 1: 创建 QTL 对象
qtl_global <- makeqtl(cross = my_cross,
                      chr = 1,
                      pos = 190.9419,
                      what = "prob")

# Step 2: 拟合模型
fit_result <- fitqtl(cross = my_cross,
                     qtl = qtl_global,
                     pheno.col = "global_methy_lvl",
                     method = "hk")  # 或者 "imp", 看你有没有做 imputation

# Step 3: 查看结果
summary(fit_result)


scan_global[250:350,]
post <- fread("/mnt/int/RIL/for_paper/codes_github/data/newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt")[,1:4]

start_global <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2543",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_global <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2690",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_global$lod)



scan_globalm <- scan_global[scan_global$lod==max(scan_global$lod),]
scan_globalm$annotate_label <- paste0("lod:",round(scan_globalm$lod,2),"\nInt:",start_global,"-",end_global,"Mb")
scan_globalm 


lodint(scan_gbM, chr=1, drop=2)
# chr      pos      lod      category
# 1_2623      1 184.9254 11.02450 gbM_methy_lvl -> 1_2623
# c1.loc189   1 189.0000 13.11037 gbM_methy_lvl -> 1_2655 
# 1_2673      1 193.9272 11.00076 gbM_methy_lvl ->1_2673   


scan_gbM[250:350,]
start_gbM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2623",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_gbM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2673",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_gbM$lod)



scan_gbMm <- scan_gbM[scan_gbM$lod==max(scan_gbM$lod),]
scan_gbMm$annotate_label <- paste0("lod:",round(scan_gbMm$lod,2),"\nInt:",start_gbM,"-",end_gbM,"Mb")
scan_gbMm 



library(qtl)

# Step 1: 创建 QTL 对象
qtl_gbM <- makeqtl(cross = my_cross,
                   chr = 1,
                   pos = 189.0000,
                   what = "prob")

# Step 2: 拟合模型
fit_result_gbM <- fitqtl(cross = my_cross,
                         qtl = qtl_gbM,
                         pheno.col = "gbM_methy_lvl",
                         method = "hk")  # 或 "imp" 视你的数据而定

# Step 3: 输出结果，查看 %var（即 x%）
summary(fit_result_gbM)







lodint(scan_gene, chr=1, drop=2)
# chr      pos      lod       category
# c1.loc185   1 185.0000 10.75444 gene_methy_lvl ->1_2623
# c1.loc190   1 190.0000 12.86411 gene_methy_lvl
# 1_2673      1 193.9272 10.81485 gene_methy_lvl



scan_gene[250:350,]
start_gene <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2623",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_gene <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2673",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_gene$lod)



scan_genem <- scan_gene[scan_gene$lod==max(scan_gene$lod),]
scan_genem$annotate_label <- paste0("lod:",round(scan_genem$lod,2),"\nInt:",start_gene,"-",end_gene,"Mb")
scan_genem 


# 创建 gene 表型的主 QTL 对象
qtl_gene <- makeqtl(cross = my_cross,
                    chr = 1,
                    pos = 190.000,
                    what = "prob")

# 拟合模型
fit_result_gene <- fitqtl(cross = my_cross,
                          qtl = qtl_gene,
                          pheno.col = "gene_methy_lvl",
                          method = "hk")

# 输出解释度
summary(fit_result_gene)





lodint(scan_UM, chr=1, drop=2)
# chr      pos      lod     category
# 1_2601   1 183.4327 5.691666 UM_methy_lvl
# 1_2661   1 190.9419 7.908228 UM_methy_lvl
# 1_2690   1 196.9125 5.711966 UM_methy_lvl
scan_UM[250:350,]
start_UM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2601",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_UM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="1_2690",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_UM$lod)

library(qtl)

# Step 1: 创建 QTL 对象
qtl_UM <- makeqtl(cross = my_cross,
                  chr = 1,
                  pos = 190.9419,
                  what = "prob")

# Step 2: 拟合模型
fit_result_UM <- fitqtl(cross = my_cross,
                        qtl = qtl_UM,
                        pheno.col = "UM_methy_lvl",
                        method = "hk")

# Step 3: 输出结果，查看 %var（解释度）
summary(fit_result_UM)


scan_UMm <- scan_UM[scan_UM$lod==max(scan_UM$lod),]
scan_UMm$annotate_label <- paste0("lod:",round(scan_UMm$lod,2),"\nInt:",start_UM,"-",end_UM,"Mb")
scan_UMm 


lodint(scan_TE, chr=2, drop=2)
# chr      pos      lod     category
# 2_2267   2 122.3804 1.635201 TE_methy_lvl
# 2_2357   2 140.3269 3.875601 TE_methy_lvl 
# 2_2377   2 141.8196 3.365601 TE_methy_lvl


scan_TE[250:350,]
start_TE <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2267",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_TE <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_TE$lod)


library(qtl)

# Step 1: 创建 QTL 对象
qtl_TE <- makeqtl(cross = my_cross,
                  chr = 2,
                  pos = 140.3269,
                  what = "prob")

# Step 2: 拟合模型
fit_result_TE <- fitqtl(cross = my_cross,
                        qtl = qtl_TE,
                        pheno.col = "TE_methy_lvl",
                        method = "hk")  # 或 "imp"，看你有没有 imputation

# Step 3: 输出结果，查看 %var（x%）
summary(fit_result_TE)




scan_TEm <- scan_TE[scan_TE$lod==max(scan_TE$lod),]
scan_TEm$annotate_label <- paste0("lod:",round(scan_TEm$lod,2),"\nInt:",start_TE,"-",end_TE,"Mb")
scan_TEm 


lodint(scan_intergenic, chr=2, drop=2)
# chr      pos      lod             category
# 2_2207      2 111.2134 1.945974 intergenic_methy_lvl
# c2.loc131   2 131.0000 4.376524 intergenic_methy_lvl. ->2_2315
# 2_2377      2 141.8196 3.418372 intergenic_methy_lvl
scan_intergenic[650:850,]


library(qtl)

# Step 1: 创建 QTL 对象
qtl_intergenic <- makeqtl(cross = my_cross,
                          chr = 2,
                          pos = 131.0000,
                          what = "prob")

# Step 2: 拟合模型
fit_result_intergenic <- fitqtl(cross = my_cross,
                                qtl = qtl_intergenic,
                                pheno.col = "intergenic_methy_lvl",
                                method = "hk")  # 可替换为 "imp" 如果做过imputation

# Step 3: 查看解释度
summary(fit_result_intergenic)


start_intergenic <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2207",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_intergenic <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_intergenic$lod)



scan_intergenicm <- scan_intergenic[scan_intergenic$lod==max(scan_intergenic$lod),]
scan_intergenicm$annotate_label <- paste0("lod:",round(scan_intergenicm$lod,2),"\nInt:",start_intergenic,"-",end_intergenic,"Mb")
scan_intergenicm 






scan_teM <- scanone(my_cross, pheno.col = "teM_methy_lvl", method = "hk")
#scan_gbM <- scanone(my_cross, pheno.col = "gbM_methy_lvl", method = "hk")



scan_teM$category <- "teM_methy_lvl"
#scan_gbM$category <- "gbM_methy_lvl"



perm_results_teM <- scanone(my_cross, pheno.col = "teM_methy_lvl", n.perm = 1000)
su_perm_results_teM <-summary(perm_results_teM)
su_perm_results_teM
cutoff_teM <- su_perm_results_teM[1,1]



library(ggplot2)
library(scales)





# lodint(scan_gbM, chr=1, drop=2)

#chr2
lodint(scan_teM, chr=2, drop=2)
# chr      pos      lod      category
# c2.loc123   2 123.0000 3.157968 teM_methy_lvl ->2_2267
# 2_2357      2 140.3269 5.295889 teM_methy_lvl
# 2_2377      2 141.8196 4.488945 teM_methy_lvl



scan_teM[scan_teM$chr==2,]


start_teM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2267",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_teM <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_teM[scan_teM$chr==2,]$lod)

scan_teMc2 <-scan_teM[scan_teM$chr==2,]
scan_teMmc2 <- scan_teMc2[scan_teMc2$lod==max(scan_teM[scan_teM$chr==2,]$lod),]
scan_teMmc2$annotate_label <- paste0("lod:",round(scan_teMmc2$lod,2),"\nInt:",start_teM,"-",end_teM,"Mb")
scan_teMmc2


# Step 1: 创建 QTL 对象
qtl_teM <- makeqtl(cross = my_cross,
                   chr = 2,
                   pos = 140.3269,
                   what = "prob")

# Step 2: 拟合模型
fit_result_teM <- fitqtl(cross = my_cross,
                         qtl = qtl_teM,
                         pheno.col = "teM_methy_lvl",
                         method = "hk")  # 或 "imp"，看你有没有 imputation

# Step 3: 输出结果，查看 %var（x%）
summary(fit_result_teM)





scan_all <-rbind(scan_global,scan_gbM,scan_gene,scan_UM,scan_TE,scan_intergenic,scan_teM)

scan_all$category <- factor(scan_all$category,c("global_methy_lvl","gene_methy_lvl","gbM_methy_lvl","UM_methy_lvl","TE_methy_lvl","teM_methy_lvl","intergenic_methy_lvl"))

output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 10
plot.width <- 25
out.name <-"Fig2C_global_gbM_gene_UM_TE_intergenic teM lod curves new methy level 2026-01"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(scan_all, aes(x=pos, y=lod)) +  
  geom_line(aes(color=category))+
  scale_color_manual(values=c("black","#663300","blue","purple","pink","darkgreen","grey"))+
  labs(title="new methy level") +
  geom_hline(yintercept = cutoff_global, linetype = "dashed", color = "black", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_gbM, linetype = "dashed", color = "blue", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_gene, linetype = "dashed", color ="#663300", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_UM, linetype = "dashed", color = "purple", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_TE, linetype = "dashed", color = "pink", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_intergenic, linetype = "dashed", color = "grey", size = 0.2) +  # 添加水平线
  geom_hline(yintercept = cutoff_teM, linetype = "dashed", color ="darkgreen", size = 0.2) +  # 添加水平线
  
  geom_text(data =scan_globalm, aes(label = annotate_label), x =scan_globalm$pos-30, y = scan_globalm$lod-0.3, hjust = 1, vjust =1, color ="black") +
  geom_text(data =scan_UMm, aes(label = annotate_label), x =scan_UMm$pos-30, y = scan_UMm$lod-0.3, hjust = 1, vjust =1, color ="purple") +
  geom_text(data =scan_TEm, aes(label = annotate_label), x =scan_TEm$pos-20, y = scan_TEm$lod-0.2, hjust = 1, vjust =1, color ="pink") +
  geom_text(data =scan_intergenicm, aes(label = annotate_label), x =scan_intergenicm$pos-10, y = scan_intergenicm$lod-0.2, hjust = 1, vjust =1, color ="grey") +
  geom_text(data =scan_teMmc2, aes(label = annotate_label), x =scan_teMmc2$pos-10, y = scan_teMmc2$lod-0.2, hjust = 1, vjust =1, color ="darkgreen") +

  geom_text(data =scan_gbMm, aes(label = annotate_label), x =scan_gbMm$pos-20, y = scan_gbMm$lod-0.3, hjust = 1, vjust =1, color ="blue") +

  geom_text(data =scan_genem, aes(label = annotate_label), x =scan_genem$pos-20, y = scan_genem$lod-0.7, hjust = 1, vjust =1, color ="#663300") +
 
  theme_classic() +
  facet_wrap(~ chr, nrow=1, scales = 'free_x') +
  scale_y_continuous(expand = c(0,0)) +  # 使用 label_comma 避免科学计数法
  theme(text = element_text(size = 15),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )   # 添加次要网格线
print(gg)
dev.off()


