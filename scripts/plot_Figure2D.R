
#Fig2D
# ==============================================================================
# Script Name: plot_Figure2D.R
# Description: Generates Figure 2D - QTL LOD curve non-CG methylation level RILs
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


#methy lvl
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

data_file<-"newrun_2024-11_68_RIL_global_CHG_CWA_CHH_nonCWA_new_methy_lvl_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2024-11.csv"
my_cross <- read.cross(file = data_file, format = "csv", genotypes = c("0", "2"))
summary(my_cross)

my_cross_est <- est.map(my_cross,error.prob=1e-07)

my_cross <-replace.map(my_cross, my_cross_est)
my_cross



my_cross<- calc.genoprob(my_cross, step=1)

scan_CWA <- scanone(my_cross, pheno.col = "global_CWA_methy_lvl", method = "hk")
scan_CHH <- scanone(my_cross, pheno.col = "global_CHH_methy_lvl", method = "hk")
scan_nonCWA <- scanone(my_cross, pheno.col = "global_nonCWA_methy_lvl", method = "hk")
scan_CHG <- scanone(my_cross, pheno.col = "global_CHG_methy_lvl", method = "hk")


scan_CWA$category <- "global_CWA_methy_lvl"
scan_CHH$category <- "global_CHH_methy_lvl"
scan_nonCWA$category <- "global_nonCWA_methy_lvl"
scan_CHG$category <- "global_CHG_methy_lvl"


perm_results_CWA <- scanone(my_cross, pheno.col = "global_CWA_methy_lvl", n.perm = 1000)
su_perm_results_CWA <-summary(perm_results_CWA)
su_perm_results_CWA
cutoff_CWA <- su_perm_results_CWA[1,1]

perm_results_CHH <- scanone(my_cross, pheno.col = "global_CHH_methy_lvl", n.perm = 1000)
su_perm_results_CHH <-summary(perm_results_CHH)
su_perm_results_CHH
cutoff_CHH <- su_perm_results_CHH[1,1]

perm_results_nonCWA <- scanone(my_cross, pheno.col = "global_nonCWA_methy_lvl", n.perm = 1000)
su_perm_results_nonCWA <-summary(perm_results_nonCWA)
su_perm_results_nonCWA
cutoff_nonCWA <- su_perm_results_nonCWA[1,1]

perm_results_CHG <- scanone(my_cross, pheno.col = "global_CHG_methy_lvl", n.perm = 1000)
su_perm_results_CHG <-summary(perm_results_CHG)
su_perm_results_CHG
cutoff_CHG <- su_perm_results_CHG[1,1]
# 
scan_all <-rbind(scan_CWA,scan_CHH,scan_nonCWA,scan_CHG)




post <- fread("/mnt/int/RIL/for_paper/codes_github/data/newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt")[,1:4]




#chr2
lodint(scan_CWA, chr=2, drop=2)
# chr      pos      lod             category
# 2_2257      2 119.3951 1.421325 global_CWA_methy_lvl
# c2.loc138   2 138.0000 3.462464 global_CWA_methy_lvl
# 2_2377      2 141.8196 2.699015 global_CWA_methy_lvl

#check major QTL was identified on chromosome 1, explaining  x% of total CG methylation level variance
library(qtl)

# Step 1: 创建 QTL 对象
qtl_CWA<- makeqtl(cross = my_cross,
                  chr = 2,
                  pos =138.0000,
                  what = "prob")

# Step 2: 拟合模型
fit_result <- fitqtl(cross = my_cross,
                     qtl = qtl_CWA,
                     pheno.col = "global_CWA_methy_lvl",
                     method = "hk")  # 或者 "imp", 看你有没有做 imputation

# Step 3: 查看结果
summary(fit_result)




scan_CWA[scan_CWA$chr==2,][200:250,]


start_CWA <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2257",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_CWA <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_CWA[scan_CWA$chr==2,]$lod)

scan_CWAc2 <-scan_CWA[scan_CWA$chr==2,]
scan_CWAmc2 <- scan_CWAc2[scan_CWAc2$lod==max(scan_CWA[scan_CWA$chr==2,]$lod),]
scan_CWAmc2$annotate_label <- paste0("lod:",round(scan_CWAmc2$lod,2),"\nInt:",start_CWA,"-",end_CWA,"Mb")
scan_CWAmc2


lodint(scan_nonCWA, chr=2, drop=2)
# chr      pos      lod                category
# c2.loc120   2 120.0000 3.333226 global_nonCWA_methy_lvl -> 2_2257 
# c2.loc137   2 137.0000 5.470869 global_nonCWA_methy_lvl
# 2_2377      2 141.8196 3.967868 global_nonCWA_methy_lvl

qtl_nonCWA<- makeqtl(cross = my_cross,
                     chr = 2,
                     pos =137.0000,
                     what = "prob")

# Step 2: 拟合模型
fit_result <- fitqtl(cross = my_cross,
                     qtl = qtl_nonCWA,
                     pheno.col = "global_nonCWA_methy_lvl",
                     method = "hk")  # 或者 "imp", 看你有没有做 imputation

# Step 3: 查看结果
summary(fit_result)


scan_nonCWA[scan_nonCWA$chr==2,][200:250,]


start_nonCWA <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2257",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_nonCWA <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_nonCWA[scan_nonCWA$chr==2,]$lod)

scan_nonCWAc2 <-scan_nonCWA[scan_nonCWA$chr==2,]
scan_nonCWAmc2 <- scan_nonCWAc2[scan_nonCWAc2$lod==max(scan_nonCWA[scan_nonCWA$chr==2,]$lod),]
scan_nonCWAmc2$annotate_label <- paste0("lod:",round(scan_nonCWAmc2$lod,2),"\nInt:",start_nonCWA,"-",end_nonCWA,"Mb")
scan_nonCWAmc2

lodint(scan_CHH, chr=2, drop=2)


qtl_CHH<- makeqtl(cross = my_cross,
                  chr = 2,
                  pos =138.0000,
                  what = "prob")

# Step 2: 拟合模型
fit_result <- fitqtl(cross = my_cross,
                     qtl = qtl_CHH,
                     pheno.col = "global_CHH_methy_lvl",
                     method = "hk")  # 或者 "imp", 看你有没有做 imputation

# Step 3: 查看结果
summary(fit_result)


scan_CHH[600:700,]




start_CHH <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2257",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_CHH <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_CHH[scan_CHH$chr==2,]$lod)

round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2353",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
#18.1 

scan_CHHm <- scan_CHH[scan_CHH$lod==max(scan_CHH[scan_CHH$chr==2,]$lod),]
scan_CHHm$annotate_label <- paste0("lod:",round(scan_CHHm$lod,2),"\nInt:",start_CHH,"-",end_CHH,"Mb")
scan_CHHm 



#chr2
lodint(scan_CHG, chr=2, drop=2)
# chr      pos       lod         category
# 2_2206   2 110.4726 0.8143273 CHG_methy_lvl
# 2_2353   2 137.2957 2.8884709 CHG_methy_lvl
# 2_2377   2 141.8196 1.8829744 CHG_methy_lvl


qtl_CHG<- makeqtl(cross = my_cross,
                  chr = 2,
                  pos =137.2957,
                  what = "prob")

# Step 2: 拟合模型
fit_result <- fitqtl(cross = my_cross,
                     qtl = qtl_CHG,
                     pheno.col = "global_CHG_methy_lvl",
                     method = "hk")  # 或者 "imp", 看你有没有做 imputation

# Step 3: 查看结果
summary(fit_result)

scan_CHG[scan_CHG$chr==2,]


start_CHG <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2206",]$tag2,"_",2)[,1],"-",2)[,2])/1000000,1)
end_CHG <- round(as.numeric(str_split_fixed(str_split_fixed(post[post$tag=="2_2377",]$tag2,"_",2)[,2],"-",2)[,2])/1000000,1)

max(scan_CHG[scan_CHG$chr==2,]$lod)

scan_CHGc2 <-scan_CHG[scan_CHG$chr==2,]
scan_CHGmc2 <- scan_CHGc2[scan_CHGc2$lod==max(scan_CHG[scan_CHG$chr==2,]$lod),]
scan_CHGmc2$annotate_label <- paste0("lod:",round(scan_CHGmc2$lod,2),"\nInt:",start_CHG,"-",end_CHG,"Mb")
scan_CHGmc2

# 

scan_all$category <- factor(scan_all$category,c("global_CHG_methy_lvl","global_CHH_methy_lvl","global_CWA_methy_lvl","global_nonCWA_methy_lvl"))

output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 10
plot.width <- 25
out.name <-"Fig2D_global_CHG_CHH_CWA_nonCWA_lod curves new methy level 2025-06"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

gg <- ggplot(scan_all, aes(x = pos, y = lod, color = category)) +  
  geom_line() +
  labs(title = "global CHG & CHH & CWA & nonCWA new methy level") +
  
  # 横线使用手动指定颜色（确保与你主图颜色一致）
  geom_hline(yintercept = cutoff_CHG, linetype = "dashed", size = 0.3, color = "#F8B62D") +    # 青绿色
  geom_hline(yintercept = cutoff_CHH, linetype = "dashed", size = 0.3, color = "#920783") +    # 深蓝紫
  geom_hline(yintercept = cutoff_CWA, linetype = "dashed", size = 0.3, color = "#E76BF3") + # 粉紫
  geom_hline(yintercept = cutoff_nonCWA, linetype = "dashed", size = 0.3, color = "#A6D854") +    # 绿色
  scale_color_manual(values = c(
    "global_CHG_methy_lvl" = "#F8B62D",       # 黄色/橙色，对应 cutoff_CWA
    "global_CHH_methy_lvl" = "#920783",       # 深紫色，对应 cutoff_CHH
    "global_CWA_methy_lvl" = "#E76BF3",    # 粉紫色，对应 cutoff_nonCWA
    "global_nonCWA_methy_lvl" = "#A6D854"        # 浅绿色，对应 cutoff_CHG
  ))+
  
  # 各自加 label（这里自动继承颜色）
  geom_text(data = scan_CWAmc2, aes(label = annotate_label, color = category),
            x = scan_CWAmc2$pos - 20, y = scan_CWAmc2$lod - 0.3, hjust = 1, vjust = 1) +
  geom_text(data = scan_nonCWAmc2, aes(label = annotate_label, color = category),
            x = scan_nonCWAmc2$pos - 20, y = scan_nonCWAmc2$lod - 0.3, hjust = 1, vjust = 1) +
  geom_text(data = scan_CHHm, aes(label = annotate_label, color = category),
            x = scan_CHHm$pos - 10, y = scan_CHHm$lod - 0.3, hjust = 1, vjust = 1) +
  geom_text(data = scan_CHGmc2, aes(label = annotate_label, color = category),
            x = scan_CHGmc2$pos - 20, y = scan_CHGmc2$lod - 0.3, hjust = 1, vjust = 1) +
  
  theme_classic() +
  facet_wrap(~ chr, nrow = 1, scales = 'free_x') +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    text = element_text(size = 15),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

print(gg)
dev.off()



