
#Fig2A
# ==============================================================================
# Script Name: plot_Figure2A.R
# Description: Generates Figure 2A - global methylation level RILs
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

globals_mean <- globals %>% group_by(line) %>% summarise(mean_methy_lvl=mean(100*Observed_methy_level))

globals_se <- globals %>%
  group_by(line) %>%
  summarise(
    mean_methy_lvl_SE = sd(100*Observed_methy_level, na.rm = TRUE) / sqrt(n())
  )


k <- merge(globals_mean,globals_se)
head(k)


k$line <- gsub("RIL", "", k$line)


asa<-k[order(k$mean_methy_lvl),]
asa$line<-factor(asa$line,levels = asa$line)

asa$up <- asa$mean_methy_lvl+asa$mean_methy_lvl_SE
asa$low <- asa$mean_methy_lvl-asa$mean_methy_lvl_SE


# 添加分组信息
asa$group <- ifelse(asa$line == "Cvi", "Cvi",
                    ifelse(asa$line == "Ler", "Ler", "RIL"))

# 标注最大、最小、Cvi、Ler
asa$label_flag <- with(asa, 
                       line %in% c("Cvi", "Ler") | 
                         mean_methy_lvl == max(mean_methy_lvl) | 
                         mean_methy_lvl == min(mean_methy_lvl))

# 构造标签内容
asa$label_text <- with(asa, ifelse(label_flag,
                                   ifelse(line %in% c("Cvi", "Ler"),
                                          paste0(line, "-0\n", round(mean_methy_lvl, 1)),
                                          paste0("RIL", line, "\n", round(mean_methy_lvl, 1))),
                                   NA))

# 控制文字水平偏移
asa$label_x_nudge <- 0
asa$label_x_nudge[asa$line == "Cvi"]  <- 1.2
asa$label_x_nudge[asa$line == "15"]   <- 2.2
asa$label_x_nudge[asa$line == "101"]  <- -2.5

# 设置文字颜色：最大/最小标红，其他按 group
asa$label_color <- with(asa, ifelse(mean_methy_lvl == max(mean_methy_lvl) | 
                                      mean_methy_lvl == min(mean_methy_lvl),
                                    "red", group))

# 设置 bar 填充颜色：最大/最小标红，其他按 group
asa$bar_fill <- with(asa, ifelse(mean_methy_lvl == max(mean_methy_lvl) | 
                                   mean_methy_lvl == min(mean_methy_lvl),
                                 "red", group))



asa$bar_height <- asa$mean_methy_lvl - 14
asa$low_rel <- asa$low - 14
asa$up_rel <- asa$up - 14
asa$label_y <- asa$bar_height + 0.05


# 输出图文件设置
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width <- 5
out.name <- "Fig2A mean_new CG_methy_lvl RIL all lines 70 include Cvi Ler CG global dot plot 2025-05"

pdf(paste0(output.dir, out.name, ".pdf"), width = plot.width, height = plot.height)

# 绘图

gg <- ggplot(asa, aes(x = line, y = bar_height, fill = bar_fill)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = low_rel, ymax = up_rel), width = 0.2, size = 0.3,alpha = 0.2) +
  geom_text(data = subset(asa, label_flag),
            aes(y = label_y, label = label_text, color = label_color),
            size = 3, vjust = -0.5, show.legend = FALSE,
            position = position_nudge(x = subset(asa, label_flag)$label_x_nudge)) +
  labs(y = "new CG methylation level (%)", x = "RIL line", 
       title = "68 RIL lines + Cvi + Ler CG global") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, 8),
    breaks = seq(0, 8, 2),
    labels = seq(14, 22, 2)
  ) +
  scale_fill_manual(values = c("Cvi" = "#00A29A",
                               "Ler" = "#1D2088",
                               #"RIL" = "black",
                               "red" = "red")) +
  scale_color_manual(values = c("Cvi" = "#00A29A",
                                "Ler" = "#1D2088",
                                #"RIL" = "black",
                                "red" = "red"))




# 改red 为D55E00

print(gg)
dev.off()
