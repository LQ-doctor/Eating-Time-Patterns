# 韦恩图：差异蛋白与所有显著相关模块蛋白的交集
# 准备差异蛋白和显著模块蛋白的列表
diff_proteins <- sig_protein_ids  # 来自差异分析的显著蛋白
sig_modules_proteins <- significant_module_proteins  # 来自所有显著模块的蛋白

# 计算交集
intersection <- intersect(diff_proteins, sig_modules_proteins)
cat("差异蛋白质与显著模块蛋白质的交集包含", length(intersection), "个蛋白质\n")
# 保存交集蛋白质列表
intersection_df <- data.frame(ProteinID = intersection)
write.csv(intersection_df, "pro_Results/Figure2C_intersection_proteins.csv", row.names = FALSE)

# 使用ggplot2和ggvenn创建韦恩图
library(ggplot2)
library(ggvenn)  # 如果没有安装，请先运行: install.packages("ggvenn")
library(eoffice)  # 导出为PPT: install.packages("eoffice")

# 准备数据
venn_data <- list(
  "DEGs" = diff_proteins,
  "Module Proteins" = sig_modules_proteins
)

# 自定义配色方案
venn_colors <- c("#4472C4", "#70AD47")

# 创建韦恩图
venn_plot <- ggvenn(
  venn_data,
  fill_color = venn_colors,
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 4,
  show_percentage = FALSE
)

# 添加标题和主题设置
venn_plot <- venn_plot + 
  ggtitle("Overlap between DEGs and Module Proteins") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

# 保存为PDF和PNG
ggsave("pro_Results/Figure2C_Venn_Diagram_ggplot.pdf", venn_plot, width = 8, height = 6)
ggsave("pro_Results/Figure2C_Venn_Diagram_ggplot.png", venn_plot, width = 8, height = 6, dpi = 300)
