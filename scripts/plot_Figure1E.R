

#Fig1E

# ==============================================================================
# Script Name: plot_Figure1E.R
# Description: Generates Figure 1E - global methylation level  MA-accessions
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

rm(list=ls())
gc()
setwd("/mnt/int/RIL/for_paper/codes_github/data/")

genotypes <- c("Col0","Kn0","Mt0","Tsu0")

genotypes 
# 初始化一个空数据框以存储结果
all_data <- data.frame()

# 循环遍历每个基因型
for (genotype in genotypes) {
  # 从Rdata文件中加载数据
  col <- dget(paste0("ABneutral_Boot_CG_estimates_", genotype, "_global_2000.Rdata"))$boot.base
  
  # 创建新列
  #col$pm <- col$PrMMinf
  col$pm <- col$PrMMinf + 0.5 * col$PrUMinf
  #col$pm3 <- 1 - col$PrUUinf
  
  # 添加基因型信息作为一个新列
  col$Genotypes <- genotype
  
  # 将数据绑定到一起
  all_data <- rbind(all_data, col)
}

# 查看结果
print(all_data)
all_datas <- all_data[,21:22]
all_datas

ma<- fread("/mnt/int/RIL/for_paper/codes_github/data/new_observed_methy_level_using_proportion_all_samples_4_MA-accession_2025-05.txt")
ma

ma[ma$genotype=="Col0",]$observed_methy_lvl
all_datas[all_datas$Genotypes=="Col0",]$pm



# 提取所有基因型
genotypes <- unique(all_datas$Genotypes)
genotypes
# 初始化结果列表
results <- data.frame(Genotype = character(), p_value = numeric(), stringsAsFactors = FALSE)

# 循环每个基因型做 t.test
for (g in genotypes) {
  #g <- "Col0"
  obs <- ma[ma$genotype == g, ]$observed_methy_lvl
  ref <- all_datas[all_datas$Genotypes == g, ]$pm
  
  # 只有当参考值存在，且obs长度大于1时才做t检验
  if (length(ref) == 1 && length(obs) > 1) {
    test <- t.test(obs, mu = ref)
    pval <- test$p.value
  } else {
    pval <- NA
  }
  
  results <- rbind(results, data.frame(Genotype = g, p_value = pval))
}

# 输出结果
print(results)


library(dplyr)

# 标准误差函数
se <- function(x) sd(x) / sqrt(length(x))

# 计算 observed 的均值和 SE
obs_summary <- ma %>%
  group_by(genotype) %>%
  summarise(mean = mean(observed_methy_lvl) * 100,
            se = se(observed_methy_lvl) * 100) %>%
  rename(Genotypes = genotype)


# predicted 数据（从 all_datas）
predicted_df <- all_datas %>%
  rename(Genotypes = Genotypes, mean = pm) %>%
  mutate(type = "predicted", se = NA, mean = mean * 100)

# observed 数据合并 SE
observed_df <- obs_summary %>%
  mutate(type = "observed")

# 合并绘图数据
dfout <- rbind(predicted_df, observed_df)

# 标准化 Genotype 名字为 "Col-0" 形式


dfout <- merge(dfout, results, by.x = "Genotypes", by.y = "Genotype", all.x = TRUE)
dfout$p_label <- ifelse(dfout$type == "observed", paste0("p=", signif(dfout$p_value, 2)), NA)

dfout
dfout$Genotypes <- gsub("0$", "-0", dfout$Genotypes)


dfout$Genotypes <- factor(dfout$Genotypes,levels = c("Kn-0","Mt-0","Tsu-0","Col-0"))


output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name<-"Fig1E methylation level global 4 accessions"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


gg <- ggplot(data = dfout, aes(x = Genotypes, y = mean, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7, alpha = 0.8) +
  geom_errorbar(data = dfout,
                aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.7), width = 0.3, size = 0.4) +
  geom_text(data = dfout[dfout$type == "observed", ],
            aes(label = p_label, y = mean + se + 1),
            position = position_dodge(0.7), size = 3, vjust = 0) +
  labs(y = "mCG Methylation level (%)", x = NULL, fill = "Type") +
  scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) +
  scale_fill_manual(values = c("observed" ="#FFB74D", "predicted" ="#81C784")) +  # 指定颜色
  theme_classic(base_size = 14) 
#scale_fill_manual(values = c("observed" = "#E76F51", "predicted" = "#2A9D8F"))

# 保存图像

print(gg)
dev.off()









