# 加载必要的包
library(ggplot2)
library(dplyr)
library(ggrepel)
library(DESeq2)
library(patchwork)
library(WGCNA)
library(readxl)
library(tidyr)
library(VennDiagram)
library(grid)
library(RColorBrewer)
library(tibble)

# 1.1 读取数据
protein_112 <- read_excel("新protein.xlsx")
lipid_169 <- read_excel("新lipid.xlsx")
diethabit <- read_excel("60个ID.xlsx")

# 1.2 处理蛋白质数据矩阵
protein_112 <- protein_112[,-1]
protein_112 <- t(protein_112)
colnames(protein_112) <- protein_112[1,]
protein_112 <- protein_112[-1,]
protein_112 <- as.data.frame(protein_112)
protein_112$ID <- rownames(protein_112)

# 1.3 处理脂质数据矩阵
colnames(lipid_169)[1] <- 'ID'
lipid_169 <- t(lipid_169)
colnames(lipid_169) <- lipid_169[1,]
lipid_169 <- lipid_169[-1,]
lipid_169 <- as.data.frame(lipid_169)
lipid_169$ID <- rownames(lipid_169)

# 1.4 合并数据
diethabit_prolip <- left_join(diethabit, lipid_169, by='ID')
diethabit_prolip <- left_join(diethabit_prolip, protein_112, by='ID')
diethabit_prolip <- as.data.frame(diethabit_prolip)

# 1.5 删除指定样本
ids_to_remove <- c("J721136", "J720838","J730244", "J720912", 
                   "S120649", "J770242","J731332", "J660116", "J660733")
diethabit_prolip <- diethabit_prolip[!diethabit_prolip$ID %in% ids_to_remove, ]

# ===============================
# 2. 蛋白组差异分析 - 修复版
# ===============================
# 2.1 准备蛋白质数据
z_protein_diet <- diethabit_prolip[,753:2718]

# 2.2 去除特定蛋白质
cols <- c('SPTN2','NDUS1','GOGA2','CERS3','ES8L2','RHCG','TX1B3','SEMG2')
z_protein_diet <- z_protein_diet[, setdiff(names(z_protein_diet), cols)]

# 2.3 确保所有数据都是数值型
z_protein_diet <- as.data.frame(lapply(z_protein_diet, function(x) {
  # 转换为数值，非数值将变为NA
  as.numeric(as.character(x))
}))

# 检查是否存在NA值并填充
if(any(is.na(z_protein_diet))) {
  cat("警告: 数据中存在NA值，使用列均值填充\n")
  # 使用列均值填充NA
  z_protein_diet <- as.data.frame(apply(z_protein_diet, 2, function(x) {
    if(any(is.na(x))) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
    return(x)
  }))
}

# 2.3 转换为矩阵并转置
mat_protein <- as.matrix(z_protein_diet)
mat_protein <- t(mat_protein)  # 现在是蛋白 × 样本的矩阵

# 确保矩阵为正整数（DESeq2要求）
# 首先将负值设为0
mat_protein[mat_protein < 0] <- 0
# 然后将所有值乘以一个因子并四舍五入，以保持相对比例
scaling_factor <- 100
mat_protein <- round(mat_protein * scaling_factor)

# 2.4 添加列名
colnames(mat_protein) <- paste0("sample_", 1:ncol(mat_protein))

# 2.5 创建分组信息数据框
diet_group <- data.frame(
  ID = colnames(mat_protein),
  group = factor(diethabit_prolip$diethabit, 
                 levels = c(3, 2, 1),
                 labels = c("Extend_EW", "Late_TRE", "Early_TRE")),
  row.names = colnames(mat_protein)
)

# 2.6 DESeq2分析
# 创建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(
  countData = round(mat_protein), 
  colData = diet_group,
  design = ~ group
)

# 运行DESeq2的标准化和差异分析
dds <- DESeq(dds)

# 提取各组间的差异分析结果
# Late vs Extend_EW
res_late_extend <- results(dds, contrast = c("group", "Late_TRE", "Extend_EW"))

# Early vs Extend_EW
res_early_extend <- results(dds, contrast = c("group", "Early_TRE", "Extend_EW"))

# Early vs Late
res_early_late <- results(dds, contrast = c("group", "Early_TRE", "Late_TRE"))

# 将结果转换为数据框
df_late_extend <- as.data.frame(res_late_extend) %>%
  rownames_to_column("proteinID") %>%
  mutate(comparison = "Late_TRE_vs_Extend_EW")

df_early_extend <- as.data.frame(res_early_extend) %>%
  rownames_to_column("proteinID") %>%
  mutate(comparison = "Early_TRE_vs_Extend_EW")

df_early_late <- as.data.frame(res_early_late) %>%
  rownames_to_column("proteinID") %>%
  mutate(comparison = "Early_TRE_vs_Late_TRE")

# 合并所有结果
all_results <- bind_rows(df_late_extend, df_early_extend, df_early_late)

# 添加显著性标记
all_results <- all_results %>%
  mutate(
    sig = ifelse(padj < 0.05, #& abs(log2FoldChange) >= log2(1.5), "Significant", "Not Significant"),
                 label = ifelse(pvalue < 0.05, "P value < 0.05", "P value >= 0.05")
    )
    
    all_results <- all_results %>%
      mutate(
        sig = ifelse(padj < 0.05 & abs(log2FoldChange) >= log2(1.5), "Significant", "Not Significant"),
        label = ifelse(pvalue < 0.05, "P value < 0.05", "P value >= 0.05")
      )
    
    # 2.7 获取每个比较组中最显著的蛋白质
    top_proteins <- all_results %>%
      filter(!is.na(padj)) %>%
      group_by(comparison) %>%
      filter(padj < 0.05 & abs(log2FoldChange) >= log2(1.5)) %>%
      mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
      arrange(desc(score)) %>%
      slice_head(n = 5) %>%
      ungroup()
    
    # 2.8 添加标记列
    all_results <- all_results %>%
      mutate(protein_comp = paste(proteinID, comparison, sep = "_"))
    
    top_proteins <- top_proteins %>%
      mutate(protein_comp = paste(proteinID, comparison, sep = "_"))
    
    all_results$size <- ifelse(all_results$protein_comp %in% top_proteins$protein_comp, 2, 1)
    
    # 2.9 分离数据
    dt <- filter(all_results, size == 1)
    top_sig <- filter(all_results, size == 2)
    
    # 2.11 改进的火山图绘制
    # 计算每个比较组的log2FC范围（用于背景柱）
    dfbar <- all_results %>%
      group_by(comparison) %>%
      summarise(y = max(log2FoldChange, na.rm = TRUE))
    
    dfbar1 <- all_results %>%
      group_by(comparison) %>%
      summarise(y = min(log2FoldChange, na.rm = TRUE))
    
    # 计算数据的实际范围
    tile_height <- 0.3  # 色块高度
    text_size <- 3.5    # 文字大小
    spacing <- 1
    y_max <- max(all_results$log2FoldChange, na.rm = TRUE)
    y_min <- min(all_results$log2FoldChange, na.rm = TRUE)
    label_position <- y_min - spacing
    
    # 创建标签数据框
    dfcol <- data.frame(
      comparison = c("Late_TRE_vs_Extend_EW", 
                     "Early_TRE_vs_Extend_EW", 
                     "Early_TRE_vs_Late_TRE"),
      y = label_position,
      label = c("Late TRE vs EW", 
                "Early TRE vs EW", 
                "Early TRE vs Late TRE")
    )
    
    # 设置组颜色
    mycol <- c("Late_TRE_vs_Extend_EW" = "#4472C4",    # 蓝色
               "Early_TRE_vs_Extend_EW" = "#ED7D31",   # 橙色
               "Early_TRE_vs_Late_TRE" = "#70AD47")    # 绿色
    
    # 创建ggplot可视化
    p <- ggplot() +
      geom_col(data = dfbar,
               aes(x = comparison, y = y),
               fill = "#dcdcdc", alpha = 0.3) +
      geom_col(data = dfbar1,
               aes(x = comparison, y = y),
               fill = "#dcdcdc", alpha = 0.3) +  
      geom_jitter(data = dt,
                  aes(x = comparison, y = log2FoldChange, color = label),
                  size = 0.7,
                  width = 0.4) +
      geom_jitter(data = top_proteins,
                  aes(x = comparison, y = log2FoldChange, color = label),
                  size = 0.8,
                  width = 0.4) +
      geom_tile(data = dfcol,
                aes(x = comparison, y = y, fill = comparison),
                height = text_size * tile_height,
                color = "black",
                alpha = 0.3) +
      geom_text(data = dfcol,
                aes(x = comparison, y = y, label = label),
                angle = 0,
                size = text_size,
                fontface = "bold",
                hjust = 0.5,
                color = "black") +
      geom_text_repel(
        data = top_proteins,
        aes(x = comparison, y = log2FoldChange, label = proteinID),
        size = text_size * 1,
        force = 1.3,
        max.overlaps = Inf,
        box.padding = 0.3,
        segment.size = 0,
        segment.alpha = 0
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        axis.title = element_text(size = text_size * 3.25, color = "black", face = "bold"),
        axis.line.y = element_line(color = "black", linewidth = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1, 0),
        legend.text = element_text(size = text_size * 3)
      ) +
      scale_color_manual(
        values = c("P value < 0.05" = "#FFB6C1",
                   "P value >= 0.05" = "#87CEFA")
      ) +
      scale_fill_manual(
        values = mycol,
        guide = "none"
      ) +
      labs(
        title = "A",
        x = "",
        y = "log2(Fold Change)",
        color = "Significance"
      ) +
      scale_y_continuous(
        limits = c(y_min - 2,
                   y_max)
      )
    
    # 2.12 显示并保存火山图
    print(p)
    
    # 使用eoffice保存到PPT
    topptx(p, 
           filename = "pro_Results/Figure2A_Protein_Volcano.pptx",
           width = 6, 
           height = 5,
           append = FALSE,
           devsize = FALSE,
           units = "in")
    ggsave("pro_Results/Figure2A_Protein_Volcano.pdf", p, width = 10, height = 8)
    ggsave("pro_Results/Figure2A_Protein_Volcano.png", p, width = 10, height = 8, dpi = 300)
    
    # 2.13 保存差异分析结果
    write.csv(all_results, "pro_Results/protein_results.csv", row.names = FALSE)
    