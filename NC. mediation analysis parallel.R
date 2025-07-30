# --------------------------------
# 10. 分离上调和下调脂质代谢物的平行中介分析
# --------------------------------
library(lavaan)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)  # 明确加载dplyr并抑制冲突警告
library(knitr)
library(kableExtra)
library(semPlot)
library(patchwork)
library(writexl)  # 用于Excel输出
library(officer)  # 用于Word输出
library(flextable) # 用于创建漂亮的表格

# 确认核心Hub代谢物变量名
common_up_hub <- c("PC(18:0/18:0)", "PE(O/16:0/18:1)")
common_down_hub <- c("TAG(52:2)_FA14:0", "TAG(56:6)_FA18:1", "TAG(56:6)_FA20:5", "TAG(56:6)_FA22:5")

# 创建一个简化的名称映射
up_name_mapping <- c(
  "PC(18:0/18:0)" = "M1",
  "PE(O/16:0/18:1)" = "M2"
)

down_name_mapping <- c(
  "TAG(52:2)_FA14:0" = "M1",
  "TAG(56:6)_FA18:1" = "M2", 
  "TAG(56:6)_FA20:5" = "M3",
  "TAG(56:6)_FA22:5" = "M4"
)

# 简化的变量名
simplified_up_names <- up_name_mapping[common_up_hub]
names(simplified_up_names) <- common_up_hub

simplified_down_names <- down_name_mapping[common_down_hub]
names(simplified_down_names) <- common_down_hub

# --------------------------------
# 为上调代谢物创建分析数据集
# --------------------------------
mediation_data_up <- data.frame(
  X = diethabit_prolip$diethabit,
  Y = diethabit_prolip$HBP
)

# 添加上调Hub代谢物列，使用简化的变量名
for (i in 1:length(common_up_hub)) {
  hub <- common_up_hub[i]
  simple_name <- simplified_up_names[i]
  mediation_data_up[[simple_name]] <- diethabit_prolip[[hub]]
}

# 确保变量类型正确 - X是膳食模式，Y是HBP
mediation_data_up$X <- factor(mediation_data_up$X, 
                              levels = c(3, 2, 1), 
                              labels = c("Extend_EW", "Late_TRE", "Early_TRE"))
mediation_data_up$Y <- factor(mediation_data_up$Y, 
                              levels = c(0, 1), 
                              labels = c("control", "case"))

# 为了分析需要，将分类变量转换为数值
mediation_data_up$X_num <- as.numeric(mediation_data_up$X) - 1 # 膳食模式: Extend_EW=0, Late_TRE=1, Early_TRE=2
mediation_data_up$Y_num <- as.numeric(mediation_data_up$Y) - 1 # HBP: control=0, case=1

# 标准化中介变量，使结果更易于解释
for (mediator in c("M1", "M2")) {
  mediation_data_up[[mediator]] <- scale(mediation_data_up[[mediator]])
}

# --------------------------------
# 为下调代谢物创建分析数据集
# --------------------------------
mediation_data_down <- data.frame(
  X = diethabit_prolip$diethabit,
  Y = diethabit_prolip$HBP
)

# 添加下调Hub代谢物列，使用简化的变量名
for (i in 1:length(common_down_hub)) {
  hub <- common_down_hub[i]
  simple_name <- simplified_down_names[i]
  mediation_data_down[[simple_name]] <- diethabit_prolip[[hub]]
}

# 确保变量类型正确 - X是膳食模式，Y是HBP
mediation_data_down$X <- factor(mediation_data_down$X, 
                                levels = c(3, 2, 1), 
                                labels = c("Extend_EW", "Late_TRE", "Early_TRE"))
mediation_data_down$Y <- factor(mediation_data_down$Y, 
                                levels = c(0, 1), 
                                labels = c("control", "case"))

# 为了分析需要，将分类变量转换为数值
mediation_data_down$X_num <- as.numeric(mediation_data_down$X) - 1 # 膳食模式: Extend_EW=0, Late_TRE=1, Early_TRE=2
mediation_data_down$Y_num <- as.numeric(mediation_data_down$Y) - 1 # HBP: control=0, case=1

# 标准化中介变量，使结果更易于解释
for (mediator in c("M1", "M2", "M3", "M4")) {
  mediation_data_down[[mediator]] <- scale(mediation_data_down[[mediator]])
}

# --------------------------------
# 上调代谢物并行中介分析模型
# --------------------------------
# 创建上调代谢物并行中介模型
up_model <- '
  # 直接效应路径
  Y_num ~ c*X_num
  
  # 中介效应路径 a
  M1 ~ a1*X_num
  M2 ~ a2*X_num
  
  # 中介效应路径 b
  Y_num ~ b1*M1 + b2*M2
  
  # 各个中介变量的间接效应
  indirect1 := a1*b1
  indirect2 := a2*b2
  
  # 总间接效应
  total_indirect := a1*b1 + a2*b2
  
  # 总效应
  total := c + (a1*b1 + a2*b2)
  
  # 各中介变量的中介比例
  prop1 := (a1*b1)/total
  prop2 := (a2*b2)/total
  
  # 总中介比例
  prop_mediated := total_indirect/total
'

# 运行上调代谢物并行中介模型
fit_up <- sem(up_model, data = mediation_data_up)

# 提取上调代谢物模型结果
results_up <- parameterEstimates(fit_up)

# 提取感兴趣的效应估计
results_up_filtered <- results_up[results_up$label %in% c("a1", "a2",
                                                          "b1", "b2",
                                                          "c", 
                                                          "indirect1", "indirect2",
                                                          "total_indirect", "total", 
                                                          "prop1", "prop2", "prop_mediated"), ]

# --------------------------------
# 下调代谢物并行中介分析模型
# --------------------------------
# 创建下调代谢物并行中介模型
down_model <- '
  # 直接效应路径
  Y_num ~ c*X_num
  
  # 中介效应路径 a
  M1 ~ a1*X_num
  M2 ~ a2*X_num
  M3 ~ a3*X_num
  M4 ~ a4*X_num
  
  # 中介效应路径 b
  Y_num ~ b1*M1 + b2*M2 + b3*M3 + b4*M4
  
  # 各个中介变量的间接效应
  indirect1 := a1*b1
  indirect2 := a2*b2
  indirect3 := a3*b3
  indirect4 := a4*b4
  
  # 总间接效应
  total_indirect := a1*b1 + a2*b2 + a3*b3 + a4*b4
  
  # 总效应
  total := c + (a1*b1 + a2*b2 + a3*b3 + a4*b4)
  
  # 各中介变量的中介比例
  prop1 := (a1*b1)/total
  prop2 := (a2*b2)/total
  prop3 := (a3*b3)/total
  prop4 := (a4*b4)/total
  
  # 总中介比例
  prop_mediated := total_indirect/total
'

# 运行下调代谢物并行中介模型
fit_down <- sem(down_model, data = mediation_data_down)

# 提取下调代谢物模型结果
results_down <- parameterEstimates(fit_down)

# 提取感兴趣的效应估计
results_down_filtered <- results_down[results_down$label %in% c("a1", "a2", "a3", "a4",
                                                                "b1", "b2", "b3", "b4",
                                                                "c", 
                                                                "indirect1", "indirect2", "indirect3", "indirect4",
                                                                "total_indirect", "total", 
                                                                "prop1", "prop2", "prop3", "prop4", "prop_mediated"), ]

# --------------------------------
# 创建上调代谢物结果汇总表格
# --------------------------------
summary_table_up <- data.frame(
  Path = results_up_filtered$label,
  Estimate = results_up_filtered$est,
  SE = results_up_filtered$se,
  Z_value = results_up_filtered$z,
  P_value = results_up_filtered$pvalue,
  CI_lower = results_up_filtered$ci.lower,
  CI_upper = results_up_filtered$ci.upper,
  stringsAsFactors = FALSE
)

# 添加显著性标记
summary_table_up$Significance <- ""
summary_table_up$Significance[summary_table_up$P_value < 0.05] <- "*"
summary_table_up$Significance[summary_table_up$P_value < 0.01] <- "**"
summary_table_up$Significance[summary_table_up$P_value < 0.001] <- "***"

# 保存汇总表格
write.csv(summary_table_up, "Lipid_ML_Results/data/up_parallel_mediation_summary.csv", row.names = FALSE)

# --------------------------------
# 创建下调代谢物结果汇总表格
# --------------------------------
summary_table_down <- data.frame(
  Path = results_down_filtered$label,
  Estimate = results_down_filtered$est,
  SE = results_down_filtered$se,
  Z_value = results_down_filtered$z,
  P_value = results_down_filtered$pvalue,
  CI_lower = results_down_filtered$ci.lower,
  CI_upper = results_down_filtered$ci.upper,
  stringsAsFactors = FALSE
)

# 添加显著性标记
summary_table_down$Significance <- ""
summary_table_down$Significance[summary_table_down$P_value < 0.05] <- "*"
summary_table_down$Significance[summary_table_down$P_value < 0.01] <- "**"
summary_table_down$Significance[summary_table_down$P_value < 0.001] <- "***"

# 保存汇总表格
write.csv(summary_table_down, "Lipid_ML_Results/data/down_parallel_mediation_summary.csv", row.names = FALSE)

# --------------------------------
# 创建上调代谢物美观结果表格
# --------------------------------
formatted_table_up <- data.frame(
  Effect = character(nrow(summary_table_up)),
  Mediator = character(nrow(summary_table_up)),
  Estimate_CI = character(nrow(summary_table_up)),
  Proportion = character(nrow(summary_table_up)),
  stringsAsFactors = FALSE
)

# 填充Effect列
formatted_table_up$Effect[grep("^a", summary_table_up$Path)] <- "Path a (X → M)"
formatted_table_up$Effect[grep("^b", summary_table_up$Path)] <- "Path b (M → Y)"
formatted_table_up$Effect[summary_table_up$Path == "c"] <- "Direct Effect (X → Y)"
formatted_table_up$Effect[grep("^indirect[0-9]", summary_table_up$Path)] <- "Indirect Effect"
formatted_table_up$Effect[summary_table_up$Path == "total_indirect"] <- "Total Indirect Effect"
formatted_table_up$Effect[summary_table_up$Path == "total"] <- "Total Effect"
formatted_table_up$Effect[grep("^prop", summary_table_up$Path)] <- "Proportion Mediated"

# 填充Mediator列
for (i in 1:nrow(summary_table_up)) {
  path <- summary_table_up$Path[i]
  if (grepl("1$", path) | path == "prop1") {
    formatted_table_up$Mediator[i] <- common_up_hub[1]
  } else if (grepl("2$", path) | path == "prop2") {
    formatted_table_up$Mediator[i] <- common_up_hub[2]
  } else {
    formatted_table_up$Mediator[i] <- "Overall"
  }
}

# 填充Estimate_CI列
for (i in 1:nrow(summary_table_up)) {
  formatted_table_up$Estimate_CI[i] <- sprintf("%.3f (%.3f, %.3f)%s", 
                                               summary_table_up$Estimate[i], 
                                               summary_table_up$CI_lower[i], 
                                               summary_table_up$CI_upper[i], 
                                               summary_table_up$Significance[i])
}

# 填充Proportion列
for (i in 1:nrow(summary_table_up)) {
  if (grepl("^prop", summary_table_up$Path[i])) {
    formatted_table_up$Proportion[i] <- sprintf("%.1f%%", summary_table_up$Estimate[i] * 100)
  } else {
    formatted_table_up$Proportion[i] <- NA
  }
}

# --------------------------------
# 创建下调代谢物美观结果表格
# --------------------------------
formatted_table_down <- data.frame(
  Effect = character(nrow(summary_table_down)),
  Mediator = character(nrow(summary_table_down)),
  Estimate_CI = character(nrow(summary_table_down)),
  Proportion = character(nrow(summary_table_down)),
  stringsAsFactors = FALSE
)

# 填充Effect列
formatted_table_down$Effect[grep("^a", summary_table_down$Path)] <- "Path a (X → M)"
formatted_table_down$Effect[grep("^b", summary_table_down$Path)] <- "Path b (M → Y)"
formatted_table_down$Effect[summary_table_down$Path == "c"] <- "Direct Effect (X → Y)"
formatted_table_down$Effect[grep("^indirect[0-9]", summary_table_down$Path)] <- "Indirect Effect"
formatted_table_down$Effect[summary_table_down$Path == "total_indirect"] <- "Total Indirect Effect"
formatted_table_down$Effect[summary_table_down$Path == "total"] <- "Total Effect"
formatted_table_down$Effect[grep("^prop", summary_table_down$Path)] <- "Proportion Mediated"

# 填充Mediator列
for (i in 1:nrow(summary_table_down)) {
  path <- summary_table_down$Path[i]
  if (grepl("1$", path) | path == "prop1") {
    formatted_table_down$Mediator[i] <- common_down_hub[1]
  } else if (grepl("2$", path) | path == "prop2") {
    formatted_table_down$Mediator[i] <- common_down_hub[2]
  } else if (grepl("3$", path) | path == "prop3") {
    formatted_table_down$Mediator[i] <- common_down_hub[3]
  } else if (grepl("4$", path) | path == "prop4") {
    formatted_table_down$Mediator[i] <- common_down_hub[4]
  } else {
    formatted_table_down$Mediator[i] <- "Overall"
  }
}

# 填充Estimate_CI列
for (i in 1:nrow(summary_table_down)) {
  formatted_table_down$Estimate_CI[i] <- sprintf("%.3f (%.3f, %.3f)%s", 
                                                 summary_table_down$Estimate[i], 
                                                 summary_table_down$CI_lower[i], 
                                                 summary_table_down$CI_upper[i], 
                                                 summary_table_down$Significance[i])
}

# 填充Proportion列
for (i in 1:nrow(summary_table_down)) {
  if (grepl("^prop", summary_table_down$Path[i])) {
    formatted_table_down$Proportion[i] <- sprintf("%.1f%%", summary_table_down$Estimate[i] * 100)
  } else {
    formatted_table_down$Proportion[i] <- NA
  }
}

# --------------------------------
# 创建上调代谢物漂亮表格并保存
# --------------------------------
ft_up <- flextable(formatted_table_up)
ft_up <- set_header_labels(ft_up, 
                           Effect = "Effect Type", 
                           Mediator = "Mediator", 
                           Estimate_CI = "Estimate (95% CI)", 
                           Proportion = "Proportion of Total Effect")

# 基本样式设置
ft_up <- theme_booktabs(ft_up)
ft_up <- autofit(ft_up)
ft_up <- bold(ft_up, part = "header")
# 根据Effect列合并单元格
ft_up <- merge_v(ft_up, j = "Effect")
# 设置列宽
ft_up <- width(ft_up, width = c(2, 2.5, 3, 1.5))
# 整体样式
ft_up <- bg(ft_up, bg = "#f7f7f7", part = "header")
# 修改border函数调用方式
ft_up <- border(ft_up, part = "all", border.color = "gray")

# 保存为Word文档
doc_up <- read_docx()
doc_up <- body_add_par(doc_up, "Up-regulated Lipids Parallel Mediation Analysis Results", style = "heading 1")
doc_up <- body_add_par(doc_up, "Table 1: Effects of Dietary Pattern on HBP through Up-regulated Hub Lipid Metabolites", style = "heading 2")
doc_up <- body_add_flextable(doc_up, value = ft_up)
doc_up <- body_add_par(doc_up, "Note: * p < 0.05, ** p < 0.01, *** p < 0.001", style = "Normal")
print(doc_up, target = "Lipid_ML_Results/figures/up_parallel_mediation_analysis.docx")

# --------------------------------
# 创建下调代谢物漂亮表格并保存
# --------------------------------
ft_down <- flextable(formatted_table_down)
ft_down <- set_header_labels(ft_down, 
                             Effect = "Effect Type", 
                             Mediator = "Mediator", 
                             Estimate_CI = "Estimate (95% CI)", 
                             Proportion = "Proportion of Total Effect")

# 基本样式设置
ft_down <- theme_booktabs(ft_down)
ft_down <- autofit(ft_down)
ft_down <- bold(ft_down, part = "header")
# 根据Effect列合并单元格
ft_down <- merge_v(ft_down, j = "Effect")
# 设置列宽
ft_down <- width(ft_down, width = c(2, 2.5, 3, 1.5))
# 整体样式
ft_down <- bg(ft_down, bg = "#f7f7f7", part = "header")
# 修改border函数调用方式
ft_down <- border(ft_down, part = "all", border.color = "gray")

# 保存为Word文档
doc_down <- read_docx()
doc_down <- body_add_par(doc_down, "Down-regulated Lipids Parallel Mediation Analysis Results", style = "heading 1")
doc_down <- body_add_par(doc_down, "Table 1: Effects of Dietary Pattern on HBP through Down-regulated Hub Lipid Metabolites", style = "heading 2")
doc_down <- body_add_flextable(doc_down, value = ft_down)
doc_down <- body_add_par(doc_down, "Note: * p < 0.05, ** p < 0.01, *** p < 0.001", style = "Normal")
print(doc_down, target = "Lipid_ML_Results/figures/down_parallel_mediation_analysis.docx")

# 导出为Excel
writexl::write_xlsx(list("Up-regulated" = formatted_table_up, 
                         "Down-regulated" = formatted_table_down), 
                    "Lipid_ML_Results/figures/separated_parallel_mediation_analysis.xlsx")

# --------------------------------
# 创建上调代谢物并行中介分析路径图
# --------------------------------
# 提取各路径系数
a_paths_up <- results_up_filtered %>%
  filter(grepl("^a", label)) %>%
  arrange(label)

b_paths_up <- results_up_filtered %>%
  filter(grepl("^b", label)) %>%
  arrange(label)

c_path_up <- results_up_filtered %>%
  filter(label == "c")

indirect_effects_up <- results_up_filtered %>%
  filter(grepl("^indirect[0-9]", label)) %>%
  arrange(label)

total_indirect_up <- results_up_filtered %>%
  filter(label == "total_indirect")

total_effect_up <- results_up_filtered %>%
  filter(label == "total")

proportions_up <- results_up_filtered %>%
  filter(grepl("^prop[0-9]", label)) %>%
  arrange(label)

# 计算哪些路径是显著的
a_sig_up <- a_paths_up$pvalue < 0.05
b_sig_up <- b_paths_up$pvalue < 0.05
c_sig_up <- c_path_up$pvalue < 0.05
indirect_sig_up <- indirect_effects_up$pvalue < 0.05

# 创建节点数据 - 上调代谢物
node_df_up <- data.frame(
  name = c("Dietary\nPattern", common_up_hub, "HBP"),
  type = c("exo", rep("med", 2), "out"),
  x = c(0, rep(1, 2), 2),
  y = c(0, 0.5, -0.5, 0)
)

# 命名向量以方便后续匹配
med_names_up <- common_up_hub
names(med_names_up) <- paste0("M", 1:2)

# 创建路径数据 - 使用精确匹配避免警告
a_paths_df_up <- data.frame(
  from = rep("Dietary\nPattern", 2),
  to = med_names_up,
  coef = a_paths_up$est,
  pvalue = a_paths_up$pvalue,
  path_type = "a",
  stringsAsFactors = FALSE
)

b_paths_df_up <- data.frame(
  from = med_names_up,
  to = rep("HBP", 2),
  coef = b_paths_up$est,
  pvalue = b_paths_up$pvalue,
  path_type = "b",
  stringsAsFactors = FALSE
)

c_path_df_up <- data.frame(
  from = "Dietary\nPattern",
  to = "HBP",
  coef = c_path_up$est,
  pvalue = c_path_up$pvalue,
  path_type = "c",
  stringsAsFactors = FALSE
)

# 合并所有路径
edge_df_up <- rbind(a_paths_df_up, b_paths_df_up, c_path_df_up)

# 添加标签
edge_df_up$label <- sprintf("%.2f%s", 
                            edge_df_up$coef,
                            ifelse(edge_df_up$pvalue < 0.05, "*", ""))

# 创建更准确的节点查找函数
get_node_y <- function(node_name, node_data) {
  idx <- which(node_data$name == node_name)
  if(length(idx) == 1) {
    return(node_data$y[idx])
  } else {
    return(0) # 默认值，避免错误
  }
}

# 设置路径颜色和线条粗细
edge_df_up$color <- ifelse(edge_df_up$pvalue < 0.05, "#0072B2", "gray50")
edge_df_up$size <- ifelse(edge_df_up$pvalue < 0.05, 1.2, 0.8)

# 创建发表级别的上调代谢物并行中介分析路径图
parallel_up <- ggplot() +
  # 添加节点
  geom_rect(data = node_df_up, aes(xmin = x - 0.2, xmax = x + 0.2, 
                                   ymin = y - 0.15, ymax = y + 0.15),
            fill = "white", color = "black", size = 1) +
  # 添加节点标签
  geom_text(data = node_df_up, aes(x = x, y = y, label = name), size = 3.5) +
  # 手动为每个边创建精确的位置数据
  geom_segment(data = edge_df_up, 
               aes(
                 # 精确计算每条边的起点和终点坐标
                 x = ifelse(from == "Dietary\nPattern", 0.2, 1.2),
                 y = sapply(from, get_node_y, node_data = node_df_up),
                 xend = ifelse(to == "HBP", 1.8, 0.8),
                 yend = sapply(to, get_node_y, node_data = node_df_up),
                 color = color
               ),
               arrow = arrow(length = unit(0.2, "cm")), 
               size = edge_df_up$size) +
  # 添加路径系数标签
  geom_text(data = edge_df_up,
            aes(
              x = ifelse(from == "Dietary\nPattern", 0.5,
                         ifelse(to == "HBP", 1.5, 1)),
              y = ifelse(from == "Dietary\nPattern",
                         sapply(to, get_node_y, node_data = node_df_up) + 0.3,
                         ifelse(to == "HBP",
                                sapply(from, get_node_y, node_data = node_df_up) + 0.3,
                                sapply(from, get_node_y, node_data = node_df_up))),
              label = label,
              color = color
            ),
            size = 3.5) +
  # 添加间接效应和中介比例标签
  annotate("text", x = 1, y = 1.5, 
           label = sprintf("Total Indirect Effect = %.3f (95%% CI: %.3f, %.3f)",
                           total_indirect_up$est, 
                           total_indirect_up$ci.lower, 
                           total_indirect_up$ci.upper),
           size = 4, fontface = "bold") +
  # 添加总中介效应比例
  annotate("text", x = 1, y = 1.2, 
           label = sprintf("Mediated Effect = %.1f%% of Total Effect",
                           results_up_filtered$est[results_up_filtered$label == "prop_mediated"] * 100),
           size = 4, fontface = "bold") +
  # 设置主题
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic")) +
  # 设置绘图区域范围
  xlim(-0.5, 2.5) +
  ylim(-1.5, 2) +
  # 添加标题
  labs(
    title = "Up-regulated Lipids Parallel Mediation Analysis",
    subtitle = "Mediation Effects of Up-regulated Hub Lipid Metabolites between Dietary Pattern and HBP"
  )

# --------------------------------
# 创建下调代谢物并行中介分析路径图
# --------------------------------
# 提取各路径系数
a_paths_down <- results_down_filtered %>%
  filter(grepl("^a", label)) %>%
  arrange(label)

b_paths_down <- results_down_filtered %>%
  filter(grepl("^b", label)) %>%
  arrange(label)

c_path_down <- results_down_filtered %>%
  filter(label == "c")

indirect_effects_down <- results_down_filtered %>%
  filter(grepl("^indirect[0-9]", label)) %>%
  arrange(label)

total_indirect_down <- results_down_filtered %>%
  filter(label == "total_indirect")

total_effect_down <- results_down_filtered %>%
  filter(label == "total")

proportions_down <- results_down_filtered %>%
  filter(grepl("^prop[0-9]", label)) %>%
  arrange(label)

# 计算哪些路径是显著的
a_sig_down <- a_paths_down$pvalue < 0.05
b_sig_down <- b_paths_down$pvalue < 0.05
c_sig_down <- c_path_down$pvalue < 0.05
indirect_sig_down <- indirect_effects_down$pvalue < 0.05

# 创建节点数据 - 下调代谢物
node_df_down <- data.frame(
  name = c("Dietary\nPattern", common_down_hub, "HBP"),
  type = c("exo", rep("med", 4), "out"),
  x = c(0, rep(1, 4), 2),
  y = c(0, 1.5, 0.5, -0.5, -1.5, 0)
)

# 命名向量以方便后续匹配
med_names_down <- common_down_hub
names(med_names_down) <- paste0("M", 1:4)

# 创建路径数据 - 使用精确匹配避免警告
a_paths_df_down <- data.frame(
  from = rep("Dietary\nPattern", 4),
  to = med_names_down,
  coef = a_paths_down$est,
  pvalue = a_paths_down$pvalue,
  path_type = "a",
  stringsAsFactors = FALSE
)

b_paths_df_down <- data.frame(
  from = med_names_down,
  to = rep("HBP", 4),
  coef = b_paths_down$est,
  pvalue = b_paths_down$pvalue,
  path_type = "b",
  stringsAsFactors = FALSE
)

c_path_df_down <- data.frame(
  from = "Dietary\nPattern",
  to = "HBP",
  coef = c_path_down$est,
  pvalue = c_path_down$pvalue,
  path_type = "c",
  stringsAsFactors = FALSE
)

# 合并所有路径
edge_df_down <- rbind(a_paths_df_down, b_paths_df_down, c_path_df_down)

# 添加标签
edge_df_down$label <- sprintf("%.2f%s", 
                              edge_df_down$coef,
                              ifelse(edge_df_down$pvalue < 0.05, "*", ""))

# 设置路径颜色和线条粗细
edge_df_down$color <- ifelse(edge_df_down$pvalue < 0.05, "#0072B2", "gray50")
edge_df_down$size <- ifelse(edge_df_down$pvalue < 0.05, 1.2, 0.8)

# 创建发表级别的下调代谢物并行中介分析路径图
parallel_down <- ggplot() +
  # 添加节点
  geom_rect(data = node_df_down, aes(xmin = x - 0.2, xmax = x + 0.2, 
                                     ymin = y - 0.15, ymax = y + 0.15),
            fill = "white", color = "black", size = 1) +
  # 添加节点标签
  geom_text(data = node_df_down, aes(x = x, y = y, label = name), size = 3.5) +
  # 手动为每个边创建精确的位置数据
  geom_segment(data = edge_df_down, 
               aes(
                 # 精确计算每条边的起点和终点坐标
                 x = ifelse(from == "Dietary\nPattern", 0.2, 1.2),
                 y = sapply(from, get_node_y, node_data = node_df_down),
                 xend = ifelse(to == "HBP", 1.8, 0.8),
                 yend = sapply(to, get_node_y, node_data = node_df_down),
                 color = color
               ),
               arrow = arrow(length = unit(0.2, "cm")), 
               size = edge_df_down$size) +
  # 添加路径系数标签
  geom_text(data = edge_df_down,
            aes(
              x = ifelse(from == "Dietary\nPattern", 0.5,
                         ifelse(to == "HBP", 1.5, 1)),
              y = ifelse(from == "Dietary\nPattern",
                         sapply(to, get_node_y, node_data = node_df_down) + 0.3,
                         ifelse(to == "HBP",
                                sapply(from, get_node_y, node_data = node_df_down) + 0.3,
                                sapply(from, get_node_y, node_data = node_df_down))),
              label = label,
              color = color
            ),
            size = 3.5) +
  # 添加间接效应和中介比例标签
  annotate("text", x = 1, y = 2.5, 
           label = sprintf("Total Indirect Effect = %.3f (95%% CI: %.3f, %.3f)",
                           total_indirect_down$est, 
                           total_indirect_down$ci.lower, 
                           total_indirect_down$ci.upper),
           size = 4, fontface = "bold") +
  # 添加总中介效应比例
  annotate("text", x = 1, y = 2.2, 
           label = sprintf("Mediated Effect = %.1f%% of Total Effect",
                           results_down_filtered$est[results_down_filtered$label == "prop_mediated"] * 100),
           size = 4, fontface = "bold") +
  # 设置主题
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic")) +
  # 设置绘图区域范围
  xlim(-0.5, 2.5) +
  ylim(-2.5, 3) +
  # 添加标题
  labs(
    title = "Down-regulated Lipids Parallel Mediation Analysis",
    subtitle = "Mediation Effects of Down-regulated Hub Lipid Metabolites between Dietary Pattern and HBP"
  )

# 组合两个图形
combined_plot <- parallel_up / parallel_down +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title = "Separate Parallel Mediation Analysis for Up and Down-regulated Lipid Metabolites",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 保存组合图形
ggsave("Lipid_ML_Results/figures/combined_parallel_mediation.png", 
       combined_plot, width = 10, height = 12, dpi = 300)

