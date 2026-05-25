

# ==============================================================================
# Script Name: plot_Figure3A.R
# Description: Generates Figure 3A - gain rate loss rate for 68RIL
# Author: Zhilin Zhang
# Date: 2025-05
# ==============================================================================

#new RIL figure3 gain rate loss rate for 68RIL 2025-05


rm(list=ls())

library("stringr")
library("ggplot2")
library("dplyr")

library("data.table")

k <- fread("/mnt/int/RIL/for_paper/codes_github/data/Alpha_Beta_rate__70_RIL_lines.txt")

k$line <- gsub("RIL", "", k$line)

asa<-k[order(k$Alpha),]
asa$line<-factor(asa$line,levels = asa$line)

asa$up <- asa$Alpha+asa$Alpha_SE
asa$low <- asa$Alpha-asa$Alpha_SE
asa[asa$line=="Cvi",]$up <- 0.0007

onlyu <- asa[asa$line!="Cvi" &asa$line!="Ler",]
max(onlyu$Alpha)/min(onlyu$Alpha)
max(onlyu$Beta)/min(onlyu$Beta)

asas <- asa[asa$line!="Cvi",]
asas2 <- asa[asa$line=="Cvi",]


head(asa)

asa$colors <- "black"
asa$colors <- ifelse(asa$line %in% c( 58, 150, 182, 149), "red", asa$color)
asa$colors <- ifelse(asa$line %in% c(50,15,155,162), "purple", asa$color)

# ==== 分组与标注设置 ====
asa$Group <- "Middle Group"
asa$Group[asa$line %in% c( 58, 150, 182, 149)] <- "Lowest 4 Lines"
asa$Group[asa$line %in% c(50, 15, 155, 162)]   <- "Highest 4 Lines"
asa$Group[asa$line == "Cvi"] <- "Cvi"
asa$Group[asa$line == "Ler"] <- "Ler"
asa$Group <- factor(asa$Group, levels = c("Lowest 4 Lines", "Middle Group", "Highest 4 Lines", "Cvi", "Ler"))

# 需要标注的组
asa$label_flag <- asa$Group %in% c("Lowest 4 Lines", "Highest 4 Lines", "Cvi", "Ler")

# 找极值行（Top5/Bottom5）
min_line <- asa$line[asa$Group == "Lowest 4 Lines"][which.min(asa$Alpha[asa$Group == "Lowest 4 Lines"])]
max_line <- asa$line[asa$Group == "Highest 4 Lines"][which.max(asa$Alpha[asa$Group == "Highest 4 Lines"])]

# 设置标签内容
asa$label_text <- with(asa, ifelse(label_flag,
                                   ifelse(Group %in% c("Cvi", "Ler"),
                                          paste0(Group, "-0\n", round(Alpha / 1e-04, 1)),
                                          ifelse(line == min_line | line == max_line,
                                                 paste0("RIL", line, "\n", round(Alpha / 1e-04, 1)),
                                                 paste0("RIL", line))),
                                   NA))

# ==== 将 y 起点平移为 1 ====
asa$bar_height <- asa$Alpha / 1e-4 - 1
asa$low_rel    <- asa$low   / 1e-4 - 1
asa$up_rel     <- asa$up    / 1e-4 - 1
asa$label_y    <- asa$up_rel + 0.1

asa$label_x_nudge <- 0
asa$bar_fill      <- asa$Group
asa$label_color   <- asa$Group

# ==== 配色 ====
fill_colors <- c(
  "Lowest 4 Lines"  = "purple",  # 橘红
  "Middle Group"    = "#585959",  # 深灰
  "Highest 4 Lines" =  "red",
  "Cvi"             = "#00A29A",  # 青绿
  "Ler"             = "#1D2088"   # 深蓝
)

# ==== 输出路径与尺寸 ====
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width  <- 5
out.name <- "Fig3A_GainRate_trimmed_yaxis_2025-05"

pdf(file = paste0(output.dir, out.name, ".pdf"), width = plot.width, height = plot.height)

# ==== 绘图 ====
gg <- ggplot(asa, aes(x = line, y = bar_height, fill = bar_fill)) +
  geom_col(width = 0.7) +
  
  geom_errorbar(aes(ymin = low_rel, ymax = up_rel),
                width = 0.2, size = 0.3,  alpha = 0.2) +
  
  geom_text(data = subset(asa, label_flag),
            aes(x = line, y = label_y, label = label_text, color = label_color),
            size = 2, vjust = 0.5,
            position = position_nudge(x = subset(asa, label_flag)$label_x_nudge),
            show.legend = FALSE,
            inherit.aes = FALSE) +
  
  labs(
    y = expression(paste("Gain Rate, ", alpha, " (x", 10^-4, ")")),
    x = "RIL line",
    title = "Gain Rate in 68 RILs + Cvi + Ler"
  ) +
  
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 6.2),
    breaks = seq(0, 6, 1),
    labels = seq(1, 7, 1)  # 显示真实 gain rate 数值
  ) +
  
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = fill_colors) +
  
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

print(gg)
dev.off()





asb<-k[order(k$Beta),]

asb$line<-factor(asb$line,levels = asb$line)

asb$up <- asb$Beta+asb$Beta_SE
asb$low <- asb$Beta-asb$Beta_SE
asb[asb$line=="Cvi",]$up <- 0.0033


asbs <- asb[asb$line!="Cvi",]
asbs2 <- asb[asb$line=="Cvi",]
#asas2$up <- 0.0007


asb$colors <- "black"

# ==== 处理 Beta（Loss rate） ====
asb$Group <- "Middle Group"
asb$Group[asb$line %in% c(58, 150, 182, 149)] <- "Lowest 4 Lines"
asb$Group[asb$line %in% c(50, 15, 155, 162)]   <- "Highest 4 Lines"
asb$Group[asb$line == "Cvi"] <- "Cvi"
asb$Group[asb$line == "Ler"] <- "Ler"
asb$Group <- factor(asb$Group, levels = c("Lowest 4 Lines", "Middle Group", "Highest 4 Lines", "Cvi", "Ler"))

asb$label_flag <- asb$Group %in% c("Lowest 4 Lines", "Highest 4 Lines", "Cvi", "Ler")

min_line <- asb$line[asb$Group == "Lowest 4 Lines"][which.min(asb$Beta[asb$Group == "Lowest 4 Lines"])]
max_line <- asb$line[asb$Group == "Highest 4 Lines"][which.max(asb$Beta[asb$Group == "Highest 4 Lines"])]

asb$label_text <- with(asb, ifelse(label_flag,
                                   ifelse(Group %in% c("Cvi", "Ler"),
                                          paste0(Group, "-0\n", round(Beta / 1e-04, 1)),
                                          ifelse(line == min_line | line == max_line,
                                                 paste0("RIL", line, "\n", round(Beta / 1e-04, 1)),
                                                 paste0("RIL", line))),
                                   NA))

# ==== 平移 y 起点为 5 ====
asb$bar_height <- asb$Beta / 1e-4 - 5
asb$low_rel    <- asb$low  / 1e-4 - 5
asb$up_rel     <- asb$up   / 1e-4 - 5
asb$label_y    <- asb$up_rel + 0.1

asb$label_x_nudge <- 0
asb$bar_fill      <- asb$Group
asb$label_color   <- asb$Group

fill_colors <- c(
  "Lowest 4 Lines"  ="purple",
  "Middle Group"    = "#585959",
  "Highest 4 Lines" = "red",
  "Cvi"             = "#00A29A",
  "Ler"             = "#1D2088"
)

# ==== 输出 ====
output.dir <- "/mnt/int/RIL/for_paper/codes_github/figures/"
plot.height <- 5
plot.width  <- 5
out.name <- "Fig3A_LossRate_trimmed_yaxis_2025-05"

pdf(file = paste0(output.dir, out.name, ".pdf"), width = plot.width, height = plot.height)

gg <- ggplot(asb, aes(x = line, y = bar_height, fill = bar_fill)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = low_rel, ymax = up_rel),
                width = 0.2, size = 0.3, alpha = 0.2) +
  geom_text(data = subset(asb, label_flag),
            aes(x = line, y = label_y, label = label_text, color = label_color),
            size = 2, vjust = -0.5,
            position = position_nudge(x = subset(asb, label_flag)$label_x_nudge),
            inherit.aes = FALSE, show.legend = FALSE) +
  labs(
    y = expression(paste("Loss Rate, ", beta, " (x", 10^-4, ")")),
    x = "RIL line",
    title = "Loss Rate in 68 RILs + Cvi + Ler"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 30),           # 因为 bar 是 Beta - 5
    breaks = seq(0, 30, 5),
    labels = seq(5, 35, 5)       # 显示真实 Beta 值
  )+
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = fill_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )


print(gg)
dev.off()




