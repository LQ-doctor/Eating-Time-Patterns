# 根据 Early vs Conv 的比较结果选择差异脂质
rm(list=ls())
# 1. 读取数据
setwd("D:/TRE paper")
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


# 删除指定J730244

ids_to_remove <- c("J721136", "J720838","S360153", "J720912", "S120649", "J770242","J731332", "J660733")
diethabit_prolip <- diethabit_prolip[!diethabit_prolip$ID %in% ids_to_remove, ]

# 6. 处理脂质组数据
diethabit_prolip[,70:2718] <- apply(diethabit_prolip[,70:2718], 2, as.numeric)
z_lipid_diet <- scale(diethabit_prolip[,70:752], center=TRUE, scale=TRUE)
colnames(z_lipid_diet) <- colnames(diethabit_prolip[,70:752])



# 1. 创建表达矩阵
expr_matrix <- t(z_lipid_diet)
# 2. 创建分组信息 - 三组
group <- factor(diethabit_prolip$diethabit, 
                levels = c(3,2,1),
                labels = c("Conventional", "Late_TRE", "Early_TRE"))
table(group)   
# 3. 进行差异分析（考虑三组比较）
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 差异分析
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(
  Late_vs_Conv = Late_TRE - Conventional,
  Early_vs_Conv = Early_TRE - Conventional,
  Early_vs_Late = Early_TRE - Late_TRE,  # 这里改变了方向
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 获取所有比较的差异结果
top_late_conv <- topTable(fit2, coef = "Late_vs_Conv", number = Inf)
top_early_conv <- topTable(fit2, coef = "Early_vs_Conv", number = Inf)
top_early_late <- topTable(fit2, coef = "Early_vs_Late", number = Inf)  # 现在直接得到Early vs Late的结果

# 分别选择显著上调和下调的脂质（基于Early_vs_Late的比较）
sig_lipids <- top_early_late[top_early_late$adj.P.Val < 0.05, ]
# 现在logFC的方向已经是正确的，不需要反转
top_early_late_up <- sig_lipids[sig_lipids$logFC > 0, ]    # Early上调直接表现为正值
top_early_late_down <- sig_lipids[sig_lipids$logFC < 0, ]  # Early下调直接表现为负值

# 计算综合得分
score_up <- abs(top_early_late_up$logFC)
score_down <- abs(top_early_late_down$logFC)

# 将得分添加到原始数据框中
top_early_late_up$score <- score_up
top_early_late_down$score <- score_down

# 按照综合得分排序
top_early_conv_up_sorted <- top_early_conv_up[order(top_early_conv_up$score, decreasing = TRUE), ]
top_early_conv_down_sorted <- top_early_conv_down[order(top_early_conv_down$score, decreasing = TRUE), ]



# 选择前15个上调和前25个下调的脂质
top_15_up <- rownames(head(top_early_conv_up_sorted, 30))
top_25_down <- rownames(head(top_early_conv_down_sorted, 30))

# 合并上调和下调的脂质
top_lipids <- c(top_15_up, top_25_down)

# 提取表达矩阵
expr_top <- expr_matrix[top_lipids, ]

# 在每个组内分别进行聚类
# Conventional组聚类
conv_samples <- which(group == "Conventional")
conv_cor <- cor(expr_top[, conv_samples])
hc_conv <- hclust(as.dist(1-conv_cor), method="ward.D2")
idx_conv <- conv_samples[hc_conv$order]

# Late TRE组聚类
late_samples <- which(group == "Late_TRE")
late_cor <- cor(expr_top[, late_samples])
hc_late <- hclust(as.dist(1-late_cor), method="ward.D2")
idx_late <- late_samples[hc_late$order]

# Early TRE组聚类
early_samples <- which(group == "Early_TRE")
early_cor <- cor(expr_top[, early_samples])
hc_early <- hclust(as.dist(1-early_cor), method="ward.D2")
idx_early <- early_samples[hc_early$order]

# 合并排序后的索引
ordered_idx <- c(idx_conv, idx_late, idx_early)

# 重排样本顺序
expr_ordered <- expr_top[, ordered_idx]

# 保持原始的样本ID作为列名
colnames(expr_ordered) <- diethabit_prolip$ID_index[ordered_idx]
# 创建注释信息
annotation_col <- data.frame(
  TRE_Pattern = factor(group[ordered_idx], 
                       levels = c("Conventional", "Late_TRE", "Early_TRE")),
  row.names = colnames(expr_ordered)
)

# 设置注释颜色
annotation_colors <- list(
  TRE_Pattern = c("Conventional" = "#87CEFA",   
                  "Late_TRE" = "#FFF2AF",      
                  "Early_TRE" = "#FFB6C1")     
)

# 绘制热图
p1_TRE <- pheatmap(
  expr_ordered,
  cluster_rows = TRUE,      
  cluster_cols = FALSE,     
  scale = "row",           
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Top 60 Differential Lipids by Dietary Habit",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 10,
  fontsize_row = 6,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-4, 4, length.out = 101),
  border_color = NA,
  gaps_col = c(length(idx_conv), length(idx_conv) + length(idx_late)),
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation"
)

# 查看每组样本的具体顺序
print("\nConventional组样本顺序：")
print(idx_conv)
print("\nLate TRE组样本顺序：")
print(idx_late)
print("\nEarly TRE组样本顺序：")
print(idx_early)

topptx(filename = "PCA1.pptx",
       width = 6, # 输出图形的宽度和高度
       height = 5,append = FALSE,devsize = FALSE,
       units = "in" )# 图形宽度和高度的单位)


# PCA分析使用显著差异脂质
sig_expr_matrix <- expr_matrix[rownames(sig_lipids), ]
pca_result <- prcomp(t(sig_expr_matrix))
summ <- summary(pca_result)

# 创建数据框
df_pca <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group
)

# 提取PC1和PC2的方差贡献率
xlab <- paste0("PC1 (", round(summ$importance[2,1]*100, 1), "%)")
ylab <- paste0("PC2 (", round(summ$importance[2,2]*100, 1), "%)")

# 绘制PCA图
ggplot(df_pca, aes(x = PC1, y = PC2, color = Group)) +
  # 添加置信椭圆填充
  stat_ellipse(aes(fill = Group),
               type = "norm", 
               geom = "polygon",
               alpha = 0.2) +  # 移除color=NA，让椭圆显示边界线
  # 添加散点
  geom_point(size = 3) +
  # 设置标签
  labs(x = xlab, y = ylab, color = "Group") +
  # 设置颜色
  scale_fill_manual(values = c(
    "Conventional" = "#87CEFA",
    "Late_TRE" = "#FFF2AF", 
    "Early_TRE" = "#FFB6C1"
  )) +
  scale_color_manual(values = c(
    "Conventional" = "#87CEFA",
    "Late_TRE" = "#FFF2AF",
    "Early_TRE" = "#FFB6C1"
  )) +
  # 主题设置
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = FALSE) # 不显示填充的图例