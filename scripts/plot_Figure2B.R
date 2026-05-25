

#Fig2B G2 vs G9 plot

# only G2 vs only G9

# ==============================================================================
# Script Name: plot_Figure2B.R
# Description: Generates Figure 2B - global methylation level only G2 vs only G9 RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================


rm(list=ls())

library("stringr")
library("ggplot2")
library("dplyr")

library("data.table")



globals <- fread("/mnt/int/RIL/for_paper/codes_github/data/RIL_371_all_sample_CG_new_methy_level_without_C_M_global_2026-01.txt")
globals

globals$gen <- str_split_fixed(str_split_fixed(globals$Sample_ID,"_G",2)[,2],"_L",2)[,1] 
globals$lin <- str_split_fixed(str_split_fixed(globals$Sample_ID,"_L",2)[,2],"_",2)[,1] 
globals$sublin <- str_split_fixed(str_split_fixed(globals$Sample_ID,"_L",2)[,2],"_",2)[,2] 




library(dplyr)
library(stringr)

# ---------- STEP 1: 标准化 lin 和 sublin ----------
globals <- globals %>%
  mutate(
    lin = str_pad(lin, width = 2, pad = "0"),
    sublin = str_pad(sublin, width = 1, pad = "0"),
    group_id = paste(line, lin, sublin, sep = "_")
  )

# ---------- STEP 2: 构造 expected（sublin 原为3/4，改为1/2 以便与G9匹配） ----------
expected <- expand.grid(
  line = unique(globals$line),
  lin = c("01", "02"),
  sublin = c("3", "4"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    match_sublin = ifelse(sublin == "3", "1", "2"),
    group_id = paste(line, lin, match_sublin, sep = "_")
  )

# ---------- STEP 3: 提取 G2 和 G9 数据 ----------
g2 <- globals %>%
  filter(gen == "02", sublin %in% c("3", "4")) %>%
  mutate(
    sublin = ifelse(sublin == "3", "1", "2"),
    group_id = paste(line, lin, sublin, sep = "_")
  ) %>%
  select(group_id, G2_methy = Observed_methy_level)

g9 <- globals %>%
  filter(gen == "09", sublin %in% c("1", "2")) %>%
  select(group_id, G9_methy = Observed_methy_level)

# ---------- STEP 4: 合并 G2 和 G9 ----------
paired <- expected %>%
  left_join(g2, by = "group_id") %>%
  left_join(g9, by = "group_id")

for (i in which(is.na(paired$G2_methy))) {
  line_i <- paired$line[i]
  lin_i <- paired$lin[i]
  sublin_i <- paired$sublin[i]
  
  # Step 1: 同 lin 下另一 sublin
  candidates1 <- paired %>%
    filter(line == line_i, lin == lin_i, sublin != sublin_i, !is.na(G2_methy)) %>%
    pull(G2_methy)
  
  if (length(candidates1) > 0) {
    paired$G2_methy[i] <- candidates1[1]
  } else {
    # Step 2: 跨 lin 的相同 sublin
    alt_lin <- if (lin_i == "01") "02" else "01"
    candidates2 <- paired %>%
      filter(line == line_i, lin == alt_lin, sublin == sublin_i, !is.na(G2_methy)) %>%
      pull(G2_methy)
    
    if (length(candidates2) > 0) {
      paired$G2_methy[i] <- candidates2[1]
    } else {
      # Step 3: 同 line + gen02 中任一非 NA 值
      candidates3 <- globals %>%
        filter(line == line_i, gen == "02", !is.na(Observed_methy_level)) %>%
        pull(Observed_methy_level)
      
      if (length(candidates3) > 0) {
        paired$G2_methy[i] <- candidates3[1]
      }
    }
  }
}


for (i in which(is.na(paired$G9_methy))) {
  line_i <- paired$line[i]
  lin_i <- paired$lin[i]
  sublin_i <- paired$sublin[i]
  
  # Step 1: 同 lin 下另一 sublin
  candidates1 <- paired %>%
    filter(line == line_i, lin == lin_i, sublin != sublin_i, !is.na(G9_methy)) %>%
    pull(G9_methy)
  
  if (length(candidates1) > 0) {
    paired$G9_methy[i] <- candidates1[1]
  } else {
    # Step 2: 跨 lin 的相同 sublin
    alt_lin <- if (lin_i == "01") "02" else "01"
    candidates2 <- paired %>%
      filter(line == line_i, lin == alt_lin, sublin == sublin_i, !is.na(G9_methy)) %>%
      pull(G9_methy)
    
    if (length(candidates2) > 0) {
      paired$G9_methy[i] <- candidates2[1]
    } else {
      # Step 3: 同 line + gen09 中任一非 NA 值
      candidates3 <- globals %>%
        filter(line == line_i, gen == "09", !is.na(Observed_methy_level)) %>%
        pull(Observed_methy_level)
      
      if (length(candidates3) > 0) {
        paired$G9_methy[i] <- candidates3[1]
      }
    }
  }
}


# ---------- STEP 7: 差值计算 ----------
paired <- paired %>%
  mutate(diff = G9_methy - G2_methy)

# ---------- STEP 8: 输出 ----------
head(paired)

paired$G2_methy <- 100*paired$G2_methy
paired$G9_methy <- 100*paired$G9_methy
paired$diff <- 100*paired$diff

library(dplyr)
library(ggplot2)

# ---------- STEP 1: 过滤只保留 RIL ----------
ril_wide <- paired %>%
  filter(grepl("^RIL", line)) %>%  # 注意只保留以 RIL 开头的行
  mutate(diff = G9_methy - G2_methy)

# ---------- STEP 2: 计算统计 ----------
t_test <- t.test(ril_wide$diff, mu = 0)
mean_diff <- mean(ril_wide$diff, na.rm = TRUE)
p_value <- signif(t_test$p.value, 3)


output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name<-"Fig2B hist Comparison of mCG levels between G2 and G9 across RILs 2025-05"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

# ---------- STEP 3: 绘图 ----------
gg <- ggplot(ril_wide, aes(x =diff)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", boundary = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = mean_diff, linetype = "dotted", color = "red") +
  annotate("text", x = mean_diff, y = 56, label = paste0("mean = ", round(mean_diff, 3)),
           hjust = -0.1, color = "red", size = 4) +
  annotate("text", x = -4.8, y = 56, label = paste0("p = ", p_value),
           hjust = 0, size = 4) +
  labs(
    title = "G2 vs G9 mCG difference (RIL only)",
    x = "mCG (G9 - G2)(%)",
    y = "Count"
  ) +
  #coord_cartesian(xlim = c(-10, 10)) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))

print(gg)
dev.off()



library(dplyr)
library(tidyr)
library(ggplot2)

# 筛选 RIL
ril <- paired %>%
  filter(grepl("^RIL", line)) %>%
  rename(G02 = G2_methy, G09 = G9_methy) %>%
  mutate(diff = G09 - G02)

# 配对 t 检验
t_result <- t.test(ril$G02, ril$G09, paired = TRUE)
mean_diff <- signif(mean(ril$G02 - ril$G09), 3)
p_val <- signif(t_result$p.value, 3)

# 转换为 long 格式
ril_long <- ril %>%
  select(line, G02, G09) %>%
  pivot_longer(cols = c(G02, G09), names_to = "Generation", values_to = "mCG")

# 添加数值型 Generation 用于斜率
ril_long$Generation_num <- ifelse(ril_long$Generation == "G02", 1, 2)

# 计算均值
mean_data <- ril_long %>%
  group_by(Generation_num) %>%
  summarise(mean_mCG = mean(mCG, na.rm = TRUE)) %>%
  arrange(Generation_num)

# 斜率
slope <- signif((mean_data$mean_mCG[2] - mean_data$mean_mCG[1]), 3)


output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig2B dot plot Comparison of mCG levels between G2 and G9 across RILs 2025-05"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)

# 绘图
gg <- ggplot(ril_long, aes(x = Generation, y = mCG)) +
  geom_jitter(width = 0.15, height = 0, color = "steelblue", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", linetype = "dashed") +
  annotate("text", x = 1.5, y = max(ril_long$mCG, na.rm = TRUE)-2, 
           label = paste0("slope = ", slope), color = "red", size = 4, hjust = 0.5, vjust = -0.5) +
  labs(
    title = "mCG levels: G2 vs G9",
    subtitle = paste0("Paired t-test: mean diff = ", mean_diff, ", p = ", p_val),
    x = "Generation", y = "mCG level"
  ) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  theme_classic()

print(gg)
dev.off()





# diff lineage
rm(list=ls())

library("stringr")
library("ggplot2")
library("dplyr")

library("data.table")



globals <- fread("/mnt/int/RIL/for_paper/codes_github/data/RIL_371_all_sample_CG_new_methy_level_without_C_M_global_2026-01.txt")
globals
globals$gen <- str_split_fixed(str_split_fixed(globals$Sample_ID,"_G",2)[,2],"_L",2)[,1] 
globals$lin <- str_split_fixed(str_split_fixed(globals$Sample_ID,"_L",2)[,2],"_",2)[,1] 
globals$sublin <- str_split_fixed(str_split_fixed(globals$Sample_ID,"_L",2)[,2],"_",2)[,2] 

globals <- globals[globals$gen=="02",]
library(dplyr)
library(stringr)
library(tidyr)

# ---------- STEP 1: 标准化 lin 和 sublin ----------
globals <- globals %>%
  mutate(
    lin = str_pad(lin, width = 2, pad = "0"),
    sublin = str_pad(sublin, width = 1, pad = "0"),
    group_id = paste(line, lin, sublin, sep = "_")
  )

# ---------- STEP 2: 构造 expected（每个 RIL 2 lin × sublin 3/4 → 2 pair） ----------
expected <- expand.grid(
  line = unique(globals$line),
  sublin = c("3", "4"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    group_id_01 = paste(line, "01", sublin, sep = "_"),
    group_id_02 = paste(line, "02", sublin, sep = "_"),
    pair_id = paste(line, sublin, sep = "_")
  )

# ---------- STEP 3: 提取 G2（仅）数据 ----------
g2_data <- globals %>%
  filter(gen == "02", sublin %in% c("3", "4")) %>%
  select(group_id, Observed_methy_level)

# ---------- STEP 4: 合并到 expected ----------
paired <- expected %>%
  left_join(g2_data, by = c("group_id_01" = "group_id")) %>%
  rename(G2_01 = Observed_methy_level) %>%
  left_join(g2_data, by = c("group_id_02" = "group_id")) %>%
  rename(G2_02 = Observed_methy_level)

# ---------- STEP 5: 缺失补全（宽松第三步） ----------
for (i in 1:nrow(paired)) {
  line_i <- paired$line[i]
  sublin_i <- paired$sublin[i]
  
  # ------- 补 G2_01 -------
  if (is.na(paired$G2_01[i])) {
    # Step 1: 同 line，另一 sublin 的 G2_01
    candidates1 <- paired %>%
      filter(line == line_i, sublin != sublin_i, !is.na(G2_01)) %>%
      pull(G2_01)
    
    if (length(candidates1) > 0) {
      paired$G2_01[i] <- candidates1[1]
    } else {
      # Step 2: 同 line，当前 sublin 的 G2_02
      candidates2 <- paired %>%
        filter(line == line_i, sublin == sublin_i, !is.na(G2_02)) %>%
        pull(G2_02)
      
      if (length(candidates2) > 0) {
        paired$G2_01[i] <- candidates2[1]
      } else {
        # Step 3: globals 中所有该 line 的 G2 值，不限制 lin/sublin
        candidates3 <- globals %>%
          filter(line == line_i, gen == "02", !is.na(Observed_methy_level)) %>%
          pull(Observed_methy_level)
        
        if (length(candidates3) > 0) {
          paired$G2_01[i] <- candidates3[1]
        }
      }
    }
  }
  
  # ------- 补 G2_02 -------
  if (is.na(paired$G2_02[i])) {
    # Step 1: 同 line，另一 sublin 的 G2_02
    candidates1 <- paired %>%
      filter(line == line_i, sublin != sublin_i, !is.na(G2_02)) %>%
      pull(G2_02)
    
    if (length(candidates1) > 0) {
      paired$G2_02[i] <- candidates1[1]
    } else {
      # Step 2: 同 line，当前 sublin 的 G2_01
      candidates2 <- paired %>%
        filter(line == line_i, sublin == sublin_i, !is.na(G2_01)) %>%
        pull(G2_01)
      
      if (length(candidates2) > 0) {
        paired$G2_02[i] <- candidates2[1]
      } else {
        # Step 3: globals 中所有该 line 的 G2 值，不限制 lin/sublin
        candidates3 <- globals %>%
          filter(line == line_i, gen == "02", !is.na(Observed_methy_level)) %>%
          pull(Observed_methy_level)
        
        if (length(candidates3) > 0) {
          paired$G2_02[i] <- candidates3[1]
        }
      }
    }
  }
}

# ---------- STEP 6: 计算差值 ----------
paired <- paired %>%
  mutate(diff = G2_02 - G2_01)

# ---------- STEP 7: 输出 ----------
head(paired)

paired$G2_01 <- 100*paired$G2_01
paired$G2_02 <- 100*paired$G2_02
paired$diff <- 100*paired$diff

library(dplyr)
library(ggplot2)

# 假设 paired 是你的数据框
ril_wide <- paired %>%
  filter(grepl("^RIL", line)) %>%
  rename(G2_1 = G2_01, G2_2 = G2_02) %>%
  mutate(diff = G2_1 - G2_2)

# t 检验
t_test <- t.test(ril_wide$diff, mu = 0)
mean_diff <- mean(ril_wide$diff)
p_value <- signif(t_test$p.value, 3)

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name<-"Fig2B hist Comparison of mCG levels between G2 L1 and G2 L2 across RILs 2025-05"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 绘图
gg <- ggplot(ril_wide, aes(x = diff)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", boundary = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = mean_diff, linetype = "dotted", color = "red") +
  annotate("text", x = mean_diff, y = 25, label = paste0("mean = ", signif(mean_diff, 3)), hjust = -0.1, color = "red") +
  annotate("text", x = min(ril_wide$diff), y = 23.5, label = paste0("p = ", p_value), hjust = 0, size = 4) +
  labs(title = "G2_1 vs G2_2 mCG difference", x = "mCG (G2_1 - G2_2) (%)", y = "Count") +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  theme_classic()


print(gg)
dev.off()



library(dplyr)
library(ggplot2)
library(tidyr)



# 筛选 RIL
ril <- paired %>%
  filter(grepl("^RIL", line)) %>%
  rename(G2_1 = G2_01, G2_2 = G2_02)

# 转换为 long 格式
ril_long <- ril %>%
  select(line, G2_1, G2_2) %>%
  pivot_longer(cols = c(G2_1, G2_2), names_to = "Generation", values_to = "mCG")

# 配对 t 检验
t_test_pair <- t.test(ril$G2_1, ril$G2_2, paired = TRUE)
p_val_pair <- signif(t_test_pair$p.value, 3)
mean_diff_pair <- signif(mean(ril$G2_1 - ril$G2_2), 3)

# 数值横坐标用于斜率
ril_long$Generation_num <- ifelse(ril_long$Generation == "G2_1", 1, 2)

# 计算每组均值
mean_data <- ril_long %>%
  group_by(Generation_num) %>%
  summarise(mean_mCG = mean(mCG, na.rm = TRUE)) %>%
  arrange(Generation_num)

# 计算斜率
slope <- signif((mean_data$mean_mCG[2] - mean_data$mean_mCG[1]), 3)

output.dir<-"/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 7
out.name<-"Fig2B dot plot Comparison of mCG levels between G2 L1 and G2 L2 across RILs 2025-05"

pdf(paste(output.dir, out.name, ".pdf", sep=""), width = plot.width, height = plot.height)


# 绘图
gg <- ggplot(ril_long, aes(x = Generation, y = mCG)) +
  geom_jitter(width = 0.15, height = 0, color = "steelblue", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", linetype = "dashed") +
  geom_line(data = mean_data, aes(x = Generation_num, y = mean_mCG),
            inherit.aes = FALSE, color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 1.5, y = max(ril_long$mCG, na.rm = TRUE), 
           label = paste0("slope = ", slope), color = "red", size = 4, hjust = 0.5, vjust = -0.5) +
  labs(
    title = "mCG levels: G2_1 vs G2_2",
    subtitle = paste0("Paired t-test: mean diff = ", mean_diff_pair, ", p = ", p_val_pair),
    x = "Generation", y = "mCG level (%)"
  ) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))

print(gg)
dev.off()

