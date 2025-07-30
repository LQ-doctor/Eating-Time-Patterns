
library(dplyr)
library(lubridate)
library(reshape2)
library(patchwork)
library(MASS)
library(Matrix)
library(mvtnorm)
library(sandwich)
library(mediation)
library(sjPlot)
library(ggplot2)
library(splines)
library(quantreg)
library(ggsignif)
library(mice)
library(colorspace)
library(grid)
library(VIM)
library(lme4)
library(lmerTest)
library(agricolae)
library(cowplot)
library(lemon)
library(ggsci)
library(bruceR)
library(forestploter)
library(mgcv)
library(gratia)
library(nnet)
library(ggthemes)
library(RColorBrewer)
library(tidyverse)

library(permute)
library(lattice)
library(vegan)
library(ggrepel)
library(ggalluvial)
library(networkD3)
library(ggdensity)
library(palmerpenguins)
library(corrplot)
library(gtsummary)
library(ggridges)
library(viridisLite)
library(viridis)
library(randomForest)
library(caret)
library(skimr)
library(survival)
library(survminer)
library("scatterplot3d")
library(rgl)
library(akima)
library(plotly)
library(ggh4x)
library(readxl)
library(dplyr)
library(purrr)

library(pheatmap)
library(RColorBrewer)
library(limma)

library(linkET)
library(tidyverse)
library(linkET)
library(ggplot2)
library(reshape2)
library(cols4all)
rm(list=ls())
# 1. 读取数据

protein_112 <- read_excel("新protein.xlsx")
lipid_169 <- read_excel("新lipid.xlsx")
diethabit <- read_excel("60个ID.xlsx")
# 3. 删除第一列
protein_112 <- protein_112[,-1]
# 4. 处理蛋白质数据矩阵
protein_112 <- t(protein_112)
colnames(protein_112) <- protein_112[1,]
protein_112 <- protein_112[-1,]
protein_112 <- as.data.frame(protein_112)
protein_112$ID <- rownames(protein_112)
# 4. 处理脂质数据矩阵
colnames(lipid_169)[1] <- 'ID'
lipid_169 <- t(lipid_169)
colnames(lipid_169) <- lipid_169[1,]
lipid_169 <- lipid_169[-1,]
lipid_169 <- as.data.frame(lipid_169)
lipid_169$ID <- rownames(lipid_169)
# 5. 合并数据
diethabit_prolip <- left_join(diethabit, lipid_169, by='ID')
diethabit_prolip <- left_join(diethabit_prolip, protein_112, by='ID')
diethabit_prolip <- as.data.frame(diethabit_prolip)
# 6. 处理脂质组数据
diethabit_prolip[,70:2718] <- apply(diethabit_prolip[,70:2718], 2, as.numeric)
z_lipid_diet <- scale(diethabit_prolip[,70:752], center=TRUE, scale=TRUE)
colnames(z_lipid_diet) <- colnames(diethabit_prolip[,70:752])



# 1. 创建表达矩阵
expr_matrix <- t(z_lipid_diet)
# 2. 创建分组信息 - 三组
group <- factor(diethabit_prolip$diethabit, 
                levels = c(3, 2, 1),
                labels = c("Conventional", "Late_TRE", "Early_TRE"))

# 3. 进行差异分析（考虑三组比较）
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 差异分析
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(
  Late_vs_Conv = Late_TRE - Conventional,
  Early_vs_Conv = Early_TRE - Conventional,
  Early_vs_Late = Early_TRE- Late_TRE,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 获取所有比较的差异结果并添加行名作为ID
top_late_conv <- data.frame(
  ID = rownames(topTable(fit2, coef = "Late_vs_Conv", number = Inf)),
  topTable(fit2, coef = "Late_vs_Conv", number = Inf),
  stringsAsFactors = FALSE
)

top_early_conv <- data.frame(
  ID = rownames(topTable(fit2, coef = "Early_vs_Conv", number = Inf)),
  topTable(fit2, coef = "Early_vs_Conv", number = Inf),
  stringsAsFactors = FALSE
)

top_late_early <- data.frame(
  ID = rownames(topTable(fit2, coef = "Early_vs_Late", number = Inf)),
  topTable(fit2, coef = "Early_vs_Late", number = Inf),
  stringsAsFactors = FALSE
)

# 添加比较组信息
top_late_conv$comparison <- "Late_TRE_vs_Conventional"
top_early_conv$comparison <- "Early_TRE_vs_Conventional"
top_late_early$comparison <- "Late_TRE_vs_Early_TRE"

# 合并所有结果
results <- rbind(top_late_conv, top_early_conv, top_late_early)

# 添加上调/下调标记
results$regulation <- ifelse(results$logFC >= 0 & results$adj.P.Val < 0.05, "UP",
                             ifelse(results$logFC <= 0 & results$adj.P.Val < 0.05, "DOWN", "NS"))

# 重新排列列的顺序
results <- results %>%
  select(ID, comparison, logFC, AveExpr, t, P.Value, adj.P.Val, B, regulation)

# 将结果写入Excel文件
#write_xlsx(results, "limma_TRE_lipid_differential_analysis_results.xlsx")
# 定义脂质类别和颜色（保持不变）
# 定义脂质类别和颜色
lipidLevel <- c(
  'BMP', 'CE', 'CER', 'DAG', 'DCER', 'DGDG', 'FFA', 'HCER', 'LCER',
  'LPC', 'LPE', 'LPS', 'PC', 'PE', 'PG', 'PI', 'SM', 'TAG'
)

lipidColor <- c(
  'BMP' = '#D32626', 'CE' = '#C49C93', 'CER' = '#2277B3',
  'DAG' = '#4DB74A', 'DCER' = '#ACC5E4', 'DGDG' = '#9266AC',
  'FFA' = '#F39396', 'HCER' = '#EF7724', 'LCER' = '#D779B0',
  'LPC' = '#F4B7D0', 'LPE' = '#93CA59', 'LPS' = '#FFD700',
  'PC' = '#FFA500', 'PE' = '#FFB366', 'PG' = '#FFE5B4',
  'PI' = '#6699CC', 'SM' = '#66CC99', 'TAG' = '#B19CD9'
)

# 添加脂质类别信息
results$subclass<- case_when(
  grepl("Hex2Cer", results$ID) ~ "LCER",
  grepl("HexCer", results$ID) ~ "HCER",
  TRUE ~ toupper(gsub("[0-9].*", "", gsub("\\(.*", "", results$ID)))
)
results$subclass <- factor(results$subclass, levels = lipidLevel)

# 创建绘图函数
create_bubble_plot <- function(data, title, show_x_axis = FALSE, show_legend = FALSE) {
  # 添加显著性标记
  data$significant <- factor(data$`p-value` < 0.05, 
                             levels = c(TRUE, FALSE),
                             labels = c("P < 0.05", "P ≥ 0.05"))
  
  set.seed(666)
  p <- ggplot(data = data, 
              aes(x = subclass, y = `Log2 fold change`, 
                  size = significant,  # 修改这里
                  fill = subclass)) +
    geom_jitter(shape = 21, 
                alpha = 0.6, 
                width = 0.35,    
                height = 0) +    
    labs(x = if(show_x_axis) 'Lipid species' else '',
         y = expression(paste(log[2], '(fold change)')),
         title = title) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_fill_manual(values = lipidColor) +
    scale_size_manual(              # 修改这里
      values = c("P < 0.05" = 1.5, 
                 "P ≥ 0.05" = 1)
    ) +
    geom_hline(
      yintercept = c(0),
      color = 'grey50', linetype = 2,
      linewidth = 0.9
    ) +
    theme_test() +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      text = element_text(family = 'sans'),
      axis.title = element_text(family = 'sans', size = 14),
      axis.text = element_text(family = 'sans', colour = "black", size = 12),
      plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
      legend.box = "vertical",
      legend.position = if(show_legend) "right" else "none",
      legend.justification = "center",
      legend.box.just = "left"
    )
  if (!show_x_axis) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  } else {
    p <- p + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )
  }
  
  if (show_legend) {
    p <- p + guides(
      fill = guide_legend(order = 1, 
                          title = "Lipid species",
                          override.aes = list(size = 2),
                          ncol = 2),
      size = guide_legend(order = 2,
                          title = "Significance",
                          ncol = 1)
    )
  }
  
  return(p)
}

# 准备数据
data_late_conv <- results %>% 
  filter(comparison == "Late_TRE_vs_Conventional") %>%
  select(ID, subclass, logFC, P.Value) %>%
  rename(`Log2 fold change` = logFC, `p-value` = P.Value)

data_early_conv <- results %>%
  filter(comparison == "Early_TRE_vs_Conventional") %>%
  select(ID, subclass, logFC, P.Value) %>%
  rename(`Log2 fold change` = logFC, `p-value` = P.Value)

data_late_early <- results %>%
  filter(comparison == "Late_TRE_vs_Early_TRE") %>%
  select(ID, subclass, logFC, P.Value) %>%
  rename(`Log2 fold change` = logFC, `p-value` = P.Value)

# 创建三个图
p1 <- create_bubble_plot(data_late_conv, "Late TRE vs Conventional", FALSE, FALSE)
p2 <- create_bubble_plot(data_early_conv, "Early TRE vs Conventional", FALSE, FALSE)
p3 <- create_bubble_plot(data_late_early, "Late TRE vs Early TRE", TRUE, TRUE)
p1
p2
p3
# 使用patchwork组合图形
combined_plot <- p1 / p2 / p3
combined_plot

topptx(filename = "12.29p4气泡plot.pptx",
       width = 11.5, # 输出图形的宽度和高度
       height = 3,append = FALSE,devsize = FALSE,
       units = "in" )# 图形宽度和高度的单位)


devtools::install_github("Hy4m/linkET", force = TRUE)
library(linkET)
library(ggplot2)
library(dplyr)
library(cols4all)


# 数据预处理
diethabit_prolip[,70:2718] <- apply(diethabit_prolip[,70:2718], 2, as.numeric)
z_lipid_diet <- scale(diethabit_prolip[,70:752], center=TRUE, scale=TRUE)
colnames(z_lipid_diet) <- colnames(diethabit_prolip[,70:752])

# 1. 创建表达矩阵
expr_matrix <- t(z_lipid_diet)

# 2. 创建分组信息 - 两组 (基于HBP)
group <- factor(diethabit_prolip$HBP,
                levels = c(0, 1),
                labels = c("Normal", "Case"))

# 3. 进行差异分析
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 差异分析
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(
  Case_vs_Normal = Case - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 获取差异分析结果
top_table <- topTable(fit2, coef = "Case_vs_Normal", number = Inf)

# 分别选择上调和下调的脂质
# 上调：logFC > 0 且 显著性排序靠前
up_lipids <- top_table[top_table$logFC > 0,]
up_lipids_sorted <- up_lipids[order(up_lipids$P.Value),]
top_25_up <- rownames(head(up_lipids_sorted, 50))

# 下调：logFC < 0 且 显著性排序靠前
down_lipids <- top_table[top_table$logFC < 0,]
down_lipids_sorted <- down_lipids[order(down_lipids$P.Value),]
top_25_down <- rownames(head(down_lipids_sorted, 50))

# 合并上调和下调的脂质
top_50_lipids <- c(top_25_up, top_25_down)

# 4. 提取表达矩阵
expr_top50 <- expr_matrix[top_50_lipids,]

# 创建结果数据框
results_df <- data.frame(
  ID = rownames(top_table),
  comparison = "Case_vs_Normal",
  logFC = top_table$logFC,
  AveExpr = top_table$AveExpr,
  t = top_table$t,
  P.Value = top_table$P.Value,
  adj.P.Val = top_table$adj.P.Val,
  B = top_table$B
)

# 添加调控方向
results_df$regulation <- ifelse(results_df$logFC > 0, "up", "down")

# 根据P值排序
results_df <- results_df[order(results_df$P.Value), ]

# 筛选显著差异的结果
sig_results <- results_df[results_df$P.Value < 0.05, ]

# 保存top 50差异脂质的表达矩阵
write.csv(sig_results, file = "HBP_lipids_expression.csv")
write.csv(top_table, file = "HBP_differential_analysis_results.csv")