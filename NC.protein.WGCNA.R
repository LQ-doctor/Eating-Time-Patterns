# 使用所有蛋白质数据进行共表达分析，而不仅是差异蛋白
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

# 转换为矩阵并转置为WGCNA所需的样本×特征矩阵
all_protein_matrix <- as.matrix(z_protein_diet)
all_protein_matrix_t <- t(all_protein_matrix)  # 转置为特征×样本的矩阵

# 确保样本ID和表型数据对应
# 首先确保列名存在且与diethabit_prolip的行名对应
sample_ids <- diethabit_prolip$ID
colnames(all_protein_matrix_t) <- sample_ids

# 转置为样本×特征矩阵用于WGCNA
expr_data_all_t <- t(all_protein_matrix_t)  # 现在是样本×特征矩阵

# 调试信息
cat("样本数 (expr_data_all_t):", nrow(expr_data_all_t), "\n")
cat("特征数 (expr_data_all_t):", ncol(expr_data_all_t), "\n")

# 准备表型数据
phenotype_cols <- c("HBP", "cIMT_L", "cIMT_R", "SBP", "DBP", "R_pwv", "L_pwv", 
                    "LVDd", "LVDs", "TG", "CHOL", "HDL", "non_HDL", "LDL", "GLU")
phenotype_data <- diethabit_prolip[, phenotype_cols]

# 重要：确保表型数据与表达数据的行名一致
rownames(phenotype_data) <- diethabit_prolip$ID  # 使用样本ID作为行名
rownames(expr_data_all_t) <- diethabit_prolip$ID  # 确保表达数据也使用相同的样本ID作为行名

# 确保两个数据框的行名一致
if(!identical(rownames(expr_data_all_t), rownames(phenotype_data))) {
  cat("警告: 表达数据和表型数据的行名不一致，进行调整\n")
  # 找出共同的样本
  common_samples <- intersect(rownames(expr_data_all_t), rownames(phenotype_data))
  if(length(common_samples) > 0) {
    # 只保留共同样本
    expr_data_all_t <- expr_data_all_t[common_samples, ]
    phenotype_data <- phenotype_data[common_samples, ]
    cat("使用", length(common_samples), "个共同样本进行分析\n")
  } else {
    # 如果没有共同样本，则可能是行名格式不同导致的，尝试重新设置
    cat("未找到共同样本，尝试重新设置行名\n")
    if(nrow(expr_data_all_t) == nrow(phenotype_data)) {
      rownames(phenotype_data) <- rownames(expr_data_all_t)
      cat("行数相同，已重新设置表型数据行名\n")
    } else {
      stop("表达数据和表型数据的样本数不同，无法继续分析！")
    }
  }
}

# 确保所有表型数据都是数值型
phenotype_data <- as.data.frame(lapply(phenotype_data, function(x) {
  as.numeric(as.character(x))
}))

# 填充NA值
if(any(is.na(phenotype_data))) {
  cat("警告: 表型数据中存在NA值，使用列均值填充\n")
  phenotype_data <- as.data.frame(apply(phenotype_data, 2, function(x) {
    if(any(is.na(x))) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
    return(x)
  }))
}

# 检查是否存在缺失值
gsg <- goodSamplesGenes(expr_data_all_t, verbose = 3)
if (!gsg$allOK) {
  # 打印问题样本和基因数量
  if (sum(!gsg$goodGenes) > 0) {
    cat("移除问题基因数量:", sum(!gsg$goodGenes), "\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("移除问题样本数量:", sum(!gsg$goodSamples), "\n")
  }
  
  # 移除问题样本和基因
  expr_data_all_t <- expr_data_all_t[gsg$goodSamples, gsg$goodGenes]
  phenotype_data <- phenotype_data[rownames(expr_data_all_t), ]
}

# 3.2 选择合适的软阈值
# 查看不同软阈值的网络拓扑结构
powers <- c(1:20)
sft <- pickSoftThreshold(expr_data_all_t, powerVector = powers, verbose = 5, 
                         networkType = "unsigned")

# 绘制软阈值选择图
pdf("pro_Results/WGCNA/Figures/01_Soft_Threshold_Selection.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
# 尺度无关拓扑拟合指数与软阈值的关系
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, 
     cex = 0.9, 
     col = "red")
# 添加水平线表示R^2=0.8的阈值
abline(h = 0.8, col = "red")

# 平均连接度与软阈值的关系
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels = powers, 
     cex = 0.9, 
     col = "red")
dev.off()

# 选择软阈值
# 自动选择最小使R^2达到0.8的power值
auto_power_index <- min(which(sft$fitIndices[,2] > 0.8 & sft$fitIndices[,3] < 0))
if (length(auto_power_index) > 0 && !is.na(auto_power_index)) {
  selected_power <- sft$fitIndices[auto_power_index, 1]
} else {
  # 如果没有值能达到0.8，选择最大的R^2对应的power
  best_fit_index <- which.max(sft$fitIndices[,2])
  selected_power <- sft$fitIndices[best_fit_index, 1]
}
cat("自动选择的最佳软阈值(power):", selected_power, "\n")

# 使用以下参数，针对小样本量或特殊数据集的优化设置
selected_power <- 6  # 使用较小的软阈值
mergeCutHeight <- 0.15  # 降低模块合并阈值
minModuleSize <- 5  # 降低最小模块大小
deepSplit <- 4  # 最大灵敏度

# 3.3 构建共表达网络和模块识别
# 使用blockwiseModules函数，更高效处理大型数据集
net <- blockwiseModules(expr_data_all_t, 
                        power = selected_power,
                        TOMType = "unsigned", 
                        minModuleSize = minModuleSize, 
                        reassignThreshold = 0, 
                        mergeCutHeight = mergeCutHeight,
                        deepSplit = deepSplit,  # 增加灵敏度
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, 
                        saveTOMFileBase = "pro_Results/WGCNA/proteinTOM",
                        verbose = 3)

# 将数字标签转换为颜色标签
mergedColors <- labels2colors(net$colors)

# 检查获得了多少个模块
module_count <- length(table(mergedColors))
cat("Number of modules detected:", module_count, "\n")
print(table(mergedColors))

# 绘制树状图和模块颜色
pdf("pro_Results/WGCNA/Figures/03_Gene_Dendrogram_Module_Colors.pdf", width = 12, height = 8)
plotDendroAndColors(net$dendrograms[[1]], 
                    colors = mergedColors[net$blockGenes[[1]]],
                    groupLabels = "Module Colors", 
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()

# 获取模块的特征基因（蛋白质）
MEs <- moduleEigengenes(expr_data_all_t, 
                        colors = mergedColors, 
                        softPower = selected_power)$eigengenes

# 确保样本ID匹配
if(!identical(rownames(MEs), rownames(phenotype_data))) {
  shared_samples <- intersect(rownames(MEs), rownames(phenotype_data))
  
  if(length(shared_samples) > 0) {
    MEs <- MEs[shared_samples, ]
    phenotype_data <- phenotype_data[shared_samples, ]
  } else if(nrow(MEs) == nrow(phenotype_data)) {
    # 如果行数相同但顺序不同，尝试重新排序
    rownames(phenotype_data) <- rownames(MEs)
  } else {
    stop("无法对齐模块特征基因与表型数据的样本!")
  }
}

# 计算模块特征基因与表型的相关性
moduleTraitCor <- cor(MEs, phenotype_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

# 为热图准备文本标签
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# 为热图准备文本标签
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# 准备一个同时显示相关系数和p值显著性标记的文本矩阵
cor_pval_matrix <- matrix("", nrow=nrow(moduleTraitCor), ncol=ncol(moduleTraitCor))

for(i in 1:nrow(moduleTraitCor)) {
  for(j in 1:ncol(moduleTraitCor)) {
    # 格式化相关系数，保留两位小数
    cor_value <- sprintf("%.2f", moduleTraitCor[i,j])
    
    # 添加显著性标记
    if(moduleTraitPvalue[i,j] < 0.001) {
      cor_pval_matrix[i,j] <- paste0(cor_value, "***")
    } else if(moduleTraitPvalue[i,j] < 0.01) {
      cor_pval_matrix[i,j] <- paste0(cor_value, "**")
    } else if(moduleTraitPvalue[i,j] < 0.05) {
      cor_pval_matrix[i,j] <- paste0(cor_value, "*")
    } else {
      cor_pval_matrix[i,j] <- cor_value
    }
  }
}

# 绘制模块-表型相关性热图
pdf("pro_Results/WGCNA/Figures/05_Module_Trait_Relationships.pdf", width = 16, height = 14)
par(mar = c(8, 12, 3, 3))  # 增加左右边距，为文本腾出空间

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(phenotype_data),
  yLabels = names(MEs),  # 注意：应该使用names(MEs)而不是colnames(MEs)
  ySymbols = names(MEs), # 同上
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = cor_pval_matrix,  # 使用包含相关系数和显著性标记的矩阵
  setStdMargins = FALSE,
  cex.text = 0.65,  # 适当减小文本大小以避免重叠
  cex.lab.y = 0.8,  # Y轴标签字体大小
  xLabelsAngle = 45,  # X轴标签角度
  zlim = c(-1, 1),
  main = "Module-Trait Relationships",
  cex.main = 2,
  cex.lab = 0.8
)

# 添加图例解释显著性标记
legend("bottomright", 
       legend = c("*** p < 0.001", "** p < 0.01", "* p < 0.05"),
       cex = 0.8,
       bty = "n")

dev.off()

# 同样的修改也应用于Figure 2B
pdf("pro_Results/Figure2B_protein_modules_trait_correlations.pdf", width = 16, height = 14)
par(mar = c(8, 12, 3, 3))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(phenotype_data),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = cor_pval_matrix,
  setStdMargins = FALSE,
  cex.text = 0.65,
  cex.lab.y = 0.8,
  xLabelsAngle = 45,
  zlim = c(-1, 1),
  main = "B",
  cex.main = 2,
  cex.lab = 0.8
)

legend("bottomright", 
       legend = c("*** p < 0.001", "** p < 0.01", "* p < 0.05"),
       cex = 0.8,
       bty = "n")

dev.off()

# 调试步骤 - 检查维度
cat("expr_data_all_t维度:", dim(expr_data_all_t), "\n")
cat("MEs维度:", dim(MEs), "\n")

# 确保样本对齐
common_samples <- intersect(rownames(expr_data_all_t), rownames(MEs))
cat("共同样本数:", length(common_samples), "\n")

# 使用相同样本的数据计算模块成员关系
expr_subset <- expr_data_all_t[common_samples, ]
MEs_subset <- MEs[common_samples, ]

# 转置表达数据进行相关性计算
expr_t <- t(expr_subset)
cat("转置后expr_t维度:", dim(expr_t), "\n")
cat("MEs_subset维度:", dim(MEs_subset), "\n")
# 另一种方法：使用WGCNA中的signedKME函数
library(WGCNA)
kME <- signedKME(expr_subset, MEs_subset)
geneModuleMembership <- as.data.frame(kME)

# 计算p值
nSamples <- nrow(expr_subset)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# 3.4 提取模块成员信息
# 计算基因与模块特征基因的相关性
geneModuleMembership <- as.data.frame(cor(t(expr_data_all_t), MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(expr_data_all_t)))

# 已成功计算:
# kME <- signedKME(expr_subset, MEs_subset)
# geneModuleMembership <- as.data.frame(kME)
# nSamples <- nrow(expr_subset)
# MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# 将模块颜色与蛋白质ID关联
modules <- unique(mergedColors)
modules <- modules[modules != "grey"]  # 排除灰色模块（未分类的蛋白质）

# 将所有模块的蛋白质信息整合到一起
all_modules <- unique(mergedColors)
module_protein_map <- list()

# 创建一个数据框来存储所有模块的蛋白质
all_module_proteins_df <- data.frame(ProteinID = character(), Module = character(), stringsAsFactors = FALSE)

for(module in all_modules) {
  # 获取属于该模块的蛋白质
  module_proteins <- colnames(expr_data_all_t)[mergedColors == module]
  module_protein_map[[module]] <- module_proteins
  
  # 输出模块信息
  cat("模块", module, "包含", length(module_proteins), "个蛋白质\n")
  
  # 将这个模块的蛋白质添加到综合数据框
  if(length(module_proteins) > 0) {
    temp_df <- data.frame(
      ProteinID = module_proteins,
      Module = rep(module, length(module_proteins)),
      stringsAsFactors = FALSE
    )
    all_module_proteins_df <- rbind(all_module_proteins_df, temp_df)
  }
}

# 保存所有模块的蛋白质到一个单一的CSV文件
write.csv(all_module_proteins_df, "pro_Results/WGCNA/all_modules_proteins.csv", row.names = FALSE)

# 3.5 找出所有与表型显著相关的模块
# 为每个模块-表型对找出显著相关关系
module_trait_significant <- matrix(FALSE, nrow=nrow(moduleTraitPvalue), ncol=ncol(moduleTraitPvalue))
rownames(module_trait_significant) <- rownames(moduleTraitPvalue)
colnames(module_trait_significant) <- colnames(moduleTraitPvalue)

# 根据p值确定显著性（p < 0.05）
for(i in 1:nrow(moduleTraitPvalue)) {
  for(j in 1:ncol(moduleTraitPvalue)) {
    if(moduleTraitPvalue[i,j] < 0.05) {
      module_trait_significant[i,j] <- TRUE
    }
  }
}

# 找出至少与一个表型显著相关的模块
sig_modules <- rownames(module_trait_significant)[rowSums(module_trait_significant) > 0]
sig_modules <- sig_modules[!grepl("^MEgrey", sig_modules)]  # 排除灰色模块

cat("与至少一个表型显著相关的模块数量:", length(sig_modules), "\n")

# 保存所有显著模块与表型的显著相关关系
sig_module_trait_table <- data.frame(
  Module = rep(NA, length(sig_modules) * ncol(moduleTraitPvalue)),
  Trait = rep(NA, length(sig_modules) * ncol(moduleTraitPvalue)),
  Correlation = rep(NA, length(sig_modules) * ncol(moduleTraitPvalue)),
  Pvalue = rep(NA, length(sig_modules) * ncol(moduleTraitPvalue))
)

row_idx <- 1
for(module in sig_modules) {
  for(j in 1:ncol(moduleTraitPvalue)) {
    if(module_trait_significant[module, j]) {
      sig_module_trait_table$Module[row_idx] <- gsub("ME", "", module)
      sig_module_trait_table$Trait[row_idx] <- colnames(moduleTraitPvalue)[j]
      sig_module_trait_table$Correlation[row_idx] <- moduleTraitCor[module, j]
      sig_module_trait_table$Pvalue[row_idx] <- moduleTraitPvalue[module, j]
      row_idx <- row_idx + 1
    }
  }
}

# 移除空行
sig_module_trait_table <- sig_module_trait_table[!is.na(sig_module_trait_table$Module), ]

# 保存显著模块-表型关系
write.csv(sig_module_trait_table, "pro_Results/WGCNA/significant_module_trait_relations.csv", row.names = FALSE)

# 将所有模块的蛋白质列表保存在一个工作簿中
all_modules <- unique(mergedColors)
all_module_proteins <- list()

# 创建一个数据框来存储所有模块的蛋白质
module_protein_df <- data.frame(ProteinID = character(), Module = character(), stringsAsFactors = FALSE)

for(module in all_modules) {
  # 获取属于该模块的蛋白质
  module_proteins <- colnames(expr_data_all_t)[mergedColors == module]
  all_module_proteins[[module]] <- module_proteins
  
  # 输出模块信息
  cat("模块", module, "包含", length(module_proteins), "个蛋白质\n")
  
  # 将这个模块的蛋白质添加到数据框
  if(length(module_proteins) > 0) {
    temp_df <- data.frame(
      ProteinID = module_proteins,
      Module = rep(module, length(module_proteins)),
      stringsAsFactors = FALSE
    )
    module_protein_df <- rbind(module_protein_df, temp_df)
  }
}

# 保存所有模块的蛋白质列表到一个CSV文件
write.csv(module_protein_df, "pro_Results/WGCNA/all_module_proteins.csv", row.names = FALSE)

# 收集所有显著模块的蛋白质
significant_module_proteins <- c()
for(module in sig_modules) {
  module_color <- gsub("ME", "", module)
  module_proteins <- colnames(expr_data_all_t)[mergedColors == module_color]
  significant_module_proteins <- union(significant_module_proteins, module_proteins)
}

cat("所有显著模块包含的蛋白质总数:", length(significant_module_proteins), "\n")

# 保存所有显著模块的蛋白质列表
write.csv(data.frame(ProteinID = significant_module_proteins),
          file = "pro_Results/WGCNA/all_significant_module_proteins.csv",
          row.names = FALSE)


# 将所有模块的蛋白质信息整合到一起
all_modules <- unique(mergedColors)
module_protein_map <- list()

# 找出显著相关的模块
# 根据p值确定显著性（p < 0.05）
module_trait_significant <- matrix(FALSE, nrow=nrow(moduleTraitPvalue), ncol=ncol(moduleTraitPvalue))
rownames(module_trait_significant) <- rownames(moduleTraitPvalue)
colnames(module_trait_significant) <- colnames(moduleTraitPvalue)

for(i in 1:nrow(moduleTraitPvalue)) {
  for(j in 1:ncol(moduleTraitPvalue)) {
    if(moduleTraitPvalue[i,j] < 0.05) {
      module_trait_significant[i,j] <- TRUE
    }
  }
}

# 找出至少与一个表型显著相关的模块
sig_modules <- rownames(module_trait_significant)[rowSums(module_trait_significant) > 0]
sig_modules_colors <- gsub("ME", "", sig_modules)  # 移除ME前缀

# 创建一个数据框来存储所有模块的蛋白质
all_module_proteins_df <- data.frame(ProteinID = character(), 
                                     Module = character(), 
                                     IsSignificant = logical(),
                                     stringsAsFactors = FALSE)

for(module in all_modules) {
  # 获取属于该模块的蛋白质
  module_proteins <- colnames(expr_data_all_t)[mergedColors == module]
  module_protein_map[[module]] <- module_proteins
  
  # 判断该模块是否显著
  is_significant <- module %in% sig_modules_colors
  
  # 输出模块信息
  cat("模块", module, "包含", length(module_proteins), "个蛋白质", 
      ifelse(is_significant, "（显著）", "（非显著）"), "\n")
  
  # 将这个模块的蛋白质添加到综合数据框
  if(length(module_proteins) > 0) {
    temp_df <- data.frame(
      ProteinID = module_proteins,
      Module = rep(module, length(module_proteins)),
      IsSignificant = rep(is_significant, length(module_proteins)),
      stringsAsFactors = FALSE
    )
    all_module_proteins_df <- rbind(all_module_proteins_df, temp_df)
  }
}

