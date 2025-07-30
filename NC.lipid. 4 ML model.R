# ===============================
# 基础设置
# ===============================
# 加载必要的包
library(tidyverse)      # 数据处理和可视化
library(ggplot2)        # 绘图
library(randomForest)   # 随机森林
library(glmnet)         # LASSO和弹性网络回归
library(e1071)          # SVM
library(caret)          # 特征选择
library(pROC)           # ROC曲线分析
library(VennDiagram)    # 韦恩图
library(gridExtra)      # 组合图形
library(pheatmap)       # 热图
library(viridis)        # 高质量配色方案
library(cowplot)        # 组合图形
library(readxl)         # 读取Excel文件
library(writexl)        # 写入Excel文件
library(xgboost)        # XGBoost
library(ggpubr)         # 组合图形发布级别工具
library(scales)         # 颜色渐变
library(reshape2)       # 数据重塑
library(limma)          # 差异表达分析
rm(list=ls())
# 设置随机种子确保结果可重复
set.seed(123)

# 创建结果目录
dir.create("Lipid_ML_Results", showWarnings = FALSE)
dir.create("Lipid_ML_Results/figures", showWarnings = FALSE)
dir.create("Lipid_ML_Results/data", showWarnings = FALSE)

# ===============================
# 1. 数据读取与预处理
# ===============================
# 读取数据
protein_112 <- read_excel("新protein.xlsx")
lipid_169 <- read_excel("新lipid.xlsx")
diethabit <- read_excel("60个ID.xlsx")

# 处理蛋白质数据矩阵 
protein_112 <- protein_112[,-1]  # 删除第一列
protein_112 <- t(protein_112)
colnames(protein_112) <- protein_112[1,]
protein_112 <- protein_112[-1,]
protein_112 <- as.data.frame(protein_112)
protein_112$ID <- rownames(protein_112)

# 处理脂质数据矩阵
colnames(lipid_169)[1] <- 'ID'
lipid_169 <- t(lipid_169)
colnames(lipid_169) <- lipid_169[1,]
lipid_169 <- lipid_169[-1,]
lipid_169 <- as.data.frame(lipid_169)
lipid_169$ID <- rownames(lipid_169)

# 合并数据
diethabit_prolip <- left_join(diethabit, lipid_169, by='ID')
diethabit_prolip <- left_join(diethabit_prolip, protein_112, by='ID')
diethabit_prolip <- as.data.frame(diethabit_prolip)

# 删除指定的ID
ids_to_remove <- c("J721136", "J720838", "J730244", "J720912", "S120649", "J770242",
                   "J731332", "J660116", "J660733")
diethabit_prolip <- diethabit_prolip[!diethabit_prolip$ID %in% ids_to_remove, ]

# 处理脂质组数据 - 将字符转为数值并进行标准化
diethabit_prolip[,70:2718] <- apply(diethabit_prolip[,70:2718], 2, as.numeric)
z_lipid_diet <- scale(diethabit_prolip[,70:752], center=TRUE, scale=TRUE)
colnames(z_lipid_diet) <- colnames(diethabit_prolip[,70:752])

# 创建表达矩阵用于差异分析
expr_matrix <- t(z_lipid_diet)

# 颜色方案设置
group_colors <- c(
  "Extend_EW" = "#4472C4", # 蓝色
  "Late_TRE" = "#ED7D31",  # 橙色
  "Early_TRE" = "#70AD47"  # 绿色
)

# 组因子 - 更新组标签
group <- factor(diethabit_prolip$diethabit, 
                levels = c(3, 2, 1), 
                labels = c("Extend_EW", "Late_TRE", "Early_TRE"))

print("数据预处理完成")

# ===============================
# 2. Limma差异分析
# ===============================
# 执行差异分析
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 运行limma差异分析
fit <- lmFit(expr_matrix, design)

# 创建对比 - 更新对比名称
contrast.matrix <- makeContrasts(
  Late_vs_Extend = Late_TRE - Extend_EW,
  Early_vs_Extend = Early_TRE - Extend_EW,
  Late_vs_Early = Late_TRE - Early_TRE,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取所有比较的结果
top_late_extend <- data.frame(
  ID = rownames(topTable(fit2, coef = "Late_vs_Extend", number = Inf)),
  topTable(fit2, coef = "Late_vs_Extend", number = Inf),
  stringsAsFactors = FALSE
)

top_early_extend <- data.frame(
  ID = rownames(topTable(fit2, coef = "Early_vs_Extend", number = Inf)),
  topTable(fit2, coef = "Early_vs_Extend", number = Inf),
  stringsAsFactors = FALSE
)

top_late_early <- data.frame(
  ID = rownames(topTable(fit2, coef = "Late_vs_Early", number = Inf)),
  topTable(fit2, coef = "Late_vs_Early", number = Inf),
  stringsAsFactors = FALSE
)

# 添加比较组信息
top_late_extend$comparison <- "Late_TRE_vs_Extend_EW"
top_early_extend$comparison <- "Early_TRE_vs_Extend_EW"
top_late_early$comparison <- "Late_TRE_vs_Early_TRE"

# 合并所有结果
results <- rbind(top_late_extend, top_early_extend, top_late_early)

# 添加调控状态
results$regulation <- ifelse(results$logFC >= 0 & results$adj.P.Val < 0.05, "UP", 
                             ifelse(results$logFC <= 0 & results$adj.P.Val < 0.05, "DOWN", "NS"))

# 重新排列列的顺序
results <- results %>% 
  select(ID, comparison, logFC, AveExpr, t, P.Value, adj.P.Val, B, regulation)

print("Limma差异分析完成")

# ===============================
# 3. 筛选差异表达脂质并准备机器学习数据
# ===============================
# 从差异表达结果中选择显著差异的脂质(adj.P.Val < 0.05)
sig_lipids <- results %>% 
  filter(adj.P.Val < 0.05) %>% 
  pull(ID) %>% 
  unique()

print(paste("显著差异脂质数量:", length(sig_lipids)))

# 按照上调和下调分别筛选脂质
up_lipids <- results %>% 
  filter(adj.P.Val < 0.05, regulation == "UP") %>% 
  pull(ID) %>% 
  unique()

print(paste("上调脂质数量:", length(up_lipids)))

down_lipids <- results %>% 
  filter(adj.P.Val < 0.05, regulation == "DOWN") %>% 
  pull(ID) %>% 
  unique()

print(paste("下调脂质数量:", length(down_lipids)))

# 提取特征矩阵和目标变量
X_data <- diethabit_prolip %>%
  select(all_of(sig_lipids)) %>%
  as.data.frame()

# 分别提取上调和下调的脂质矩阵
X_up <- diethabit_prolip %>%
  select(all_of(up_lipids)) %>%
  as.data.frame()

X_down <- diethabit_prolip %>%
  select(all_of(down_lipids)) %>%
  as.data.frame()

# 设置目标变量为HBP（高血压状态）
y_data <- diethabit_prolip$HBP
y_data <- as.factor(y_data)  # 转换为因子类型

# 确保因子水平是有效的R变量名
levels(y_data) <- c("control", "case")

# 检查类别分布
print(table(y_data))

# 将特征数据与目标变量合并
all_data <- cbind(X_data, HBP = y_data)
up_data <- cbind(X_up, HBP = y_data)
down_data <- cbind(X_down, HBP = y_data)

# 查看数据维度
print(paste("全部显著差异脂质数据维度:", paste(dim(X_data), collapse=" x ")))
print(paste("上调脂质数据维度:", paste(dim(X_up), collapse=" x ")))
print(paste("下调脂质数据维度:", paste(dim(X_down), collapse=" x ")))

print("机器学习数据准备完成")

# ===============================
# 4. 分割数据集
# ===============================
# 分别创建上调和下调数据的训练集和测试集
train_index <- createDataPartition(y_data, p = 0.7, list = FALSE)

# 上调脂质
X_up_train <- X_up[train_index, ]
X_up_test <- X_up[-train_index, ]

# 下调脂质
X_down_train <- X_down[train_index, ]
X_down_test <- X_down[-train_index, ]

y_train <- y_data[train_index]
y_test <- y_data[-train_index]

# 设置交叉验证参数
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# 保存特征名称
up_feature_names <- colnames(X_up)
down_feature_names <- colnames(X_down)
print(paste("上调脂质数量:", length(up_feature_names)))
print(paste("下调脂质数量:", length(down_feature_names)))

# ===============================
# 5. 设置模型参数
# ===============================
# 设置参数以确保找到足够的特征
min_up_features_selection <- 5    # 每个模型至少选择的上调特征数
min_down_features_selection <- 20  # 每个模型至少选择的下调特征数

# ===============================
# 6. 上调脂质的特征选择分析
# ===============================
print("执行上调脂质代谢物的机器学习分析...")

# --------------------------------
# 6.1 随机森林 - 上调脂质
# --------------------------------
print("执行上调脂质的随机森林特征选择...")

# 训练随机森林模型
rf_up_model <- randomForest(
  x = X_up_train, 
  y = y_train, 
  ntree = 500, 
  importance = TRUE
)

# 提取特征重要性
rf_up_importance <- importance(rf_up_model)
rf_up_importance_df <- data.frame(
  Feature = rownames(rf_up_importance),
  Importance = rf_up_importance[, "MeanDecreaseGini"]
)
rf_up_importance_df <- rf_up_importance_df[order(rf_up_importance_df$Importance, decreasing = TRUE), ]

# 保存结果
write.csv(rf_up_importance_df, "Lipid_ML_Results/data/rf_up_features.csv", row.names = FALSE)

# 不限制前20个，而是基于重要性阈值选择特征
# 计算重要性平均值的2倍作为阈值
importance_threshold <- mean(rf_up_importance_df$Importance) * 2
rf_up_selected <- rf_up_importance_df$Feature[rf_up_importance_df$Importance > importance_threshold]

# 确保至少选择了min_up_features_selection个特征
if(length(rf_up_selected) < min_up_features_selection) {
  rf_up_selected <- head(rf_up_importance_df$Feature, min_up_features_selection)
}

print(paste("随机森林选择的上调脂质特征数量:", length(rf_up_selected)))

# --------------------------------
# 6.2 弹性网络 - 上调脂质
# --------------------------------
print("执行上调脂质的弹性网络特征选择...")

# 准备数据
x_up_train_matrix <- as.matrix(X_up_train)
y_train_numeric <- as.numeric(y_train) - 1

# 执行交叉验证
elastic_up_cv <- cv.glmnet(
  x = x_up_train_matrix,
  y = y_train_numeric,
  alpha = 0.5,  # 弹性网络参数，结合L1和L2惩罚
  family = "binomial",
  nfolds = 5
)

# 使用lambda.min作为最优参数
elastic_up_lambda <- elastic_up_cv$lambda.min

# 训练模型
elastic_up_model <- glmnet(
  x = x_up_train_matrix,
  y = y_train_numeric,
  alpha = 0.5,
  family = "binomial",
  lambda = elastic_up_lambda
)

# 提取系数
elastic_up_coef <- coef(elastic_up_model, s = elastic_up_lambda)
elastic_up_coef_df <- data.frame(
  Feature = rownames(elastic_up_coef),
  Coefficient = as.numeric(elastic_up_coef)
)
elastic_up_coef_df <- elastic_up_coef_df[-1, ]  # 移除截距
elastic_up_coef_df$AbsCoef <- abs(elastic_up_coef_df$Coefficient)
elastic_up_coef_df <- elastic_up_coef_df[order(elastic_up_coef_df$AbsCoef, decreasing = TRUE), ]

# 选择非零系数的特征
elastic_up_selected <- elastic_up_coef_df$Feature[elastic_up_coef_df$Coefficient != 0]

# 如果特征太少，选择系数最大的几个
if (length(elastic_up_selected) < min_up_features_selection) {
  elastic_up_selected <- head(elastic_up_coef_df$Feature, min_up_features_selection)
}

# 保存结果
elastic_up_importance_df <- elastic_up_coef_df[elastic_up_coef_df$Feature %in% elastic_up_selected, ]
write.csv(elastic_up_importance_df, "Lipid_ML_Results/data/elastic_up_features.csv", row.names = FALSE)

print(paste("弹性网络选择的上调脂质特征数量:", length(elastic_up_selected)))

# --------------------------------
# 6.3 SVM-RFE - 上调脂质
# --------------------------------
print("执行上调脂质的SVM-RFE特征选择...")

# 设置RFE控制参数
svmRFE_ctrl <- rfeControl(
  functions = caretFuncs,
  method = "cv",
  number = 5,
  verbose = FALSE
)

# 执行RFE
# 增加sizes参数，尝试更多的特征集合
svm_up_rfe_result <- rfe(
  x = as.data.frame(X_up_train),
  y = y_train,
  sizes = c(5, 10, 15, 20, 25, 30, 40, 50),
  rfeControl = svmRFE_ctrl,
  method = "svmRadial"
)

# 提取选择的特征
svm_up_selected <- svm_up_rfe_result$optVariables

# 如果特征太少，使用更多特征
if (length(svm_up_selected) < min_up_features_selection) {
  svm_up_selected <- svm_up_rfe_result$variables$var[1:min_up_features_selection]
}

print(paste("SVM-RFE选择的上调脂质特征数量:", length(svm_up_selected)))

# 为特征计算重要性得分
svm_up_importance <- data.frame(
  Feature = svm_up_rfe_result$variables$var,
  Importance = svm_up_rfe_result$variables$Overall
)
svm_up_importance <- svm_up_importance[order(svm_up_importance$Importance, decreasing = TRUE), ]

# 保存结果
write.csv(svm_up_importance, "Lipid_ML_Results/data/svm_rfe_up_features.csv", row.names = FALSE)

# --------------------------------
# 6.4 XGBoost - 上调脂质
# --------------------------------
print("执行上调脂质的XGBoost特征选择...")

# 准备数据
dtrain_up <- xgb.DMatrix(data = as.matrix(X_up_train), label = as.numeric(y_train) - 1)

# 设置参数
xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.1,
  max_depth = 6,
  subsample = 0.8,
  colsample_bytree = 0.8
)

# 训练模型
xgb_up_model <- xgb.train(
  params = xgb_params,
  data = dtrain_up,
  nrounds = 100,
  verbose = 0
)

# 提取特征重要性
xgb_up_importance <- xgb.importance(feature_names = colnames(X_up_train), model = xgb_up_model)

# 创建重要性数据框
if (nrow(xgb_up_importance) > 0) {
  xgb_up_importance_df <- data.frame(
    Feature = xgb_up_importance$Feature,
    Importance = xgb_up_importance$Gain
  )
  xgb_up_importance_df <- xgb_up_importance_df[order(xgb_up_importance_df$Importance, decreasing = TRUE), ]
  
  # 选择重要性高于平均值的特征
  importance_threshold <- mean(xgb_up_importance_df$Importance) * 1.5
  xgb_up_selected <- xgb_up_importance_df$Feature[xgb_up_importance_df$Importance > importance_threshold]
  
  # 如果特征太少，使用重要性前N个
  if (length(xgb_up_selected) < min_up_features_selection) {
    xgb_up_selected <- head(xgb_up_importance_df$Feature, min_up_features_selection)
  }
  
  # 保存结果
  write.csv(xgb_up_importance_df, "Lipid_ML_Results/data/newxgboost_up_features.csv", row.names = FALSE)
} else {
  # 如果XGBoost没有找到任何特征，使用随机森林的结果
  xgb_up_selected <- head(rf_up_importance_df$Feature, min_up_features_selection)
  xgb_up_importance_df <- data.frame(
    Feature = xgb_up_selected,
    Importance = seq(min_up_features_selection, 1, -1)
  )
  write.csv(xgb_up_importance_df, "Lipid_ML_Results/data/xgboost_up_featurenew.csv", row.names = FALSE)
}

print(paste("XGBoost选择的上调脂质特征数量:", length(xgb_up_selected)))
# ===============================
# 7. 下调脂质的特征选择分析
# ===============================
print("执行下调脂质代谢物的机器学习分析...")

# --------------------------------
# 7.1 随机森林 - 下调脂质
# --------------------------------
print("执行下调脂质的随机森林特征选择...")

# 训练随机森林模型
rf_down_model <- randomForest(
  x = X_down_train, 
  y = y_train, 
  ntree = 500, 
  importance = TRUE
)

# 提取特征重要性
rf_down_importance <- importance(rf_down_model)
rf_down_importance_df <- data.frame(
  Feature = rownames(rf_down_importance),
  Importance = rf_down_importance[, "MeanDecreaseGini"]
)
rf_down_importance_df <- rf_down_importance_df[order(rf_down_importance_df$Importance, decreasing = TRUE), ]

# 保存结果
write.csv(rf_down_importance_df, "Lipid_ML_Results/data/rf_down_features.csv", row.names = FALSE)

# 不限制前20个，而是基于重要性阈值选择特征
# 计算重要性平均值的2倍作为阈值
importance_threshold <- mean(rf_down_importance_df$Importance) * 2
rf_down_selected <- rf_down_importance_df$Feature[rf_down_importance_df$Importance > importance_threshold]

# 确保至少选择了min_down_features_selection个特征
if(length(rf_down_selected) < min_down_features_selection) {
  rf_down_selected <- head(rf_down_importance_df$Feature, min_down_features_selection)
}

print(paste("随机森林选择的下调脂质特征数量:", length(rf_down_selected)))

# --------------------------------
# 7.2 弹性网络 - 下调脂质
# --------------------------------
print("执行下调脂质的弹性网络特征选择...")

# 准备数据
x_down_train_matrix <- as.matrix(X_down_train)

# 执行交叉验证
elastic_down_cv <- cv.glmnet(
  x = x_down_train_matrix,
  y = y_train_numeric,
  alpha = 0.5,  # 弹性网络参数，结合L1和L2惩罚
  family = "binomial",
  nfolds = 5
)

# 使用lambda.min作为最优参数
elastic_down_lambda <- elastic_down_cv$lambda.min

# 训练模型
elastic_down_model <- glmnet(
  x = x_down_train_matrix,
  y = y_train_numeric,
  alpha = 0.5,
  family = "binomial",
  lambda = elastic_down_lambda
)

# 提取系数
elastic_down_coef <- coef(elastic_down_model, s = elastic_down_lambda)
elastic_down_coef_df <- data.frame(
  Feature = rownames(elastic_down_coef),
  Coefficient = as.numeric(elastic_down_coef)
)
elastic_down_coef_df <- elastic_down_coef_df[-1, ]  # 移除截距
elastic_down_coef_df$AbsCoef <- abs(elastic_down_coef_df$Coefficient)
elastic_down_coef_df <- elastic_down_coef_df[order(elastic_down_coef_df$AbsCoef, decreasing = TRUE), ]

# 选择非零系数的特征
elastic_down_selected <- elastic_down_coef_df$Feature[elastic_down_coef_df$Coefficient != 0]

# 如果特征太少，选择系数最大的几个
if (length(elastic_down_selected) < min_down_features_selection) {
  elastic_down_selected <- head(elastic_down_coef_df$Feature, min_down_features_selection)
}

# 保存结果
elastic_down_importance_df <- elastic_down_coef_df[elastic_down_coef_df$Feature %in% elastic_down_selected, ]
write.csv(elastic_down_importance_df, "Lipid_ML_Results/data/elastic_down_features.csv", row.names = FALSE)

print(paste("弹性网络选择的下调脂质特征数量:", length(elastic_down_selected)))

# --------------------------------
# 7.3 SVM-RFE - 下调脂质
# --------------------------------
print("执行下调脂质的SVM-RFE特征选择...")

# 执行RFE
# 增加sizes参数，尝试更多的特征集合
svm_down_rfe_result <- rfe(
  x = as.data.frame(X_down_train),
  y = y_train,
  sizes = c(5, 10, 15, 20, 25, 30, 40, 50),
  rfeControl = svmRFE_ctrl,
  method = "svmRadial"
)

# 提取选择的特征
svm_down_selected <- svm_down_rfe_result$optVariables

# 如果特征太少，使用更多特征
if (length(svm_down_selected) < min_down_features_selection) {
  # 从变量重要性列表中选择前min_down_features_selection个特征
  svm_down_selected <- svm_down_rfe_result$variables$var[1:min_down_features_selection]
  print(paste("SVM-RFE选择的特征数量不足，已扩展至", min_down_features_selection, "个特征"))
}

print(paste("SVM-RFE选择的下调脂质特征数量:", length(svm_down_selected)))

# 为特征计算重要性得分
svm_down_importance <- data.frame(
  Feature = svm_down_rfe_result$variables$var,
  Importance = svm_down_rfe_result$variables$Overall
)
svm_down_importance <- svm_down_importance[order(svm_down_importance$Importance, decreasing = TRUE), ]

# 保存结果
write.csv(svm_down_importance, "Lipid_ML_Results/data/svm_rfe_down_features.csv", row.names = FALSE)

# --------------------------------
# 7.4 XGBoost - 下调脂质
# --------------------------------
print("执行下调脂质的XGBoost特征选择...")

# 准备数据
dtrain_down <- xgb.DMatrix(data = as.matrix(X_down_train), label = as.numeric(y_train) - 1)

# 训练模型
xgb_down_model <- xgb.train(
  params = xgb_params,
  data = dtrain_down,
  nrounds = 100,
  verbose = 0
)

# 提取特征重要性
xgb_down_importance <- xgb.importance(feature_names = colnames(X_down_train), model = xgb_down_model)

# 创建重要性数据框
if (nrow(xgb_down_importance) > 0) {
  xgb_down_importance_df <- data.frame(
    Feature = xgb_down_importance$Feature,
    Importance = xgb_down_importance$Gain
  )
  xgb_down_importance_df <- xgb_down_importance_df[order(xgb_down_importance_df$Importance, decreasing = TRUE), ]
  
  # 选择重要性高于平均值的特征
  importance_threshold <- mean(xgb_down_importance_df$Importance) * 1.5
  xgb_down_selected <- xgb_down_importance_df$Feature[xgb_down_importance_df$Importance > importance_threshold]
  
  # 如果特征太少，使用重要性前N个
  if (length(xgb_down_selected) < min_down_features_selection) {
    xgb_down_selected <- head(xgb_down_importance_df$Feature, min_down_features_selection)
  }
  
  # 保存结果
  write.csv(xgb_down_importance_df, "Lipid_ML_Results/data/xgboost_down_features.csv", row.names = FALSE)
} else {
  # 如果XGBoost没有找到任何特征，使用随机森林的结果
  xgb_down_selected <- head(rf_down_importance_df$Feature, min_down_features_selection)
  xgb_down_importance_df <- data.frame(
    Feature = xgb_down_selected,
    Importance = seq(min_down_features_selection, 1, -1)
  )
  write.csv(xgb_down_importance_df, "Lipid_ML_Results/data/xgboost_down_features.csv", row.names = FALSE)
}

print(paste("XGBoost选择的下调脂质特征数量:", length(xgb_down_selected)))

# ===============================
# 8. 寻找Hub代谢物
# ===============================
print("寻找Hub代谢物...")

# --------------------------------
# 8.1 上调脂质的Hub代谢物
# --------------------------------
print("寻找上调脂质Hub代谢物...")

# 创建上调脂质各模型特征列表
up_models_list <- list(
  RF = rf_up_selected,
  Elastic = elastic_up_selected,
  SVM = svm_up_selected,
  XGBoost = xgb_up_selected
)

# 计算每个特征被选中的次数
up_all_features <- unlist(up_models_list)
up_feature_counts <- table(up_all_features)

# 选择至少被2个模型选中的特征
up_hub_features <- names(up_feature_counts[up_feature_counts >= 2])
print(paste("至少被2个模型选中的上调Hub代谢物数量:", length(up_hub_features)))

# 选择被至少3个模型选中的特征作为高置信度Hub代谢物
up_high_confidence_hub <- names(up_feature_counts[up_feature_counts >= 3])
print(paste("至少被3个模型选中的上调高置信度Hub代谢物数量:", length(up_high_confidence_hub)))

# 选择被4个模型全部选中的特征作为核心Hub代谢物
up_core_hub <- names(up_feature_counts[up_feature_counts == 4])
print(paste("被全部4个模型选中的核心上调Hub代谢物数量:", length(up_core_hub)))

# 保存上调Hub代谢物
up_hub_df <- data.frame(
  Feature = names(up_feature_counts),
  Count = as.numeric(up_feature_counts),
  Regulation = "Up-regulated"
)
up_hub_df <- up_hub_df[order(up_hub_df$Count, decreasing = TRUE), ]
write.csv(up_hub_df, "Lipid_ML_Results/data/up_hub_metabolites.csv", row.names = FALSE)

# --------------------------------
# 8.2 下调脂质的Hub代谢物
# --------------------------------
print("寻找下调脂质Hub代谢物...")

# 创建下调脂质各模型特征列表
down_models_list <- list(
  RF = rf_down_selected,
  Elastic = elastic_down_selected,
  SVM = svm_down_selected,
  XGBoost = xgb_down_selected
)

# 计算每个特征被选中的次数
down_all_features <- unlist(down_models_list)
down_feature_counts <- table(down_all_features)

# 选择至少被2个模型选中的特征
down_hub_features <- names(down_feature_counts[down_feature_counts >= 2])
print(paste("至少被2个模型选中的下调Hub代谢物数量:", length(down_hub_features)))

# 选择被至少3个模型选中的特征作为高置信度Hub代谢物
down_high_confidence_hub <- names(down_feature_counts[down_feature_counts >= 3])
print(paste("至少被3个模型选中的下调高置信度Hub代谢物数量:", length(down_high_confidence_hub)))

# 选择被4个模型全部选中的特征作为核心Hub代谢物
down_core_hub <- names(down_feature_counts[down_feature_counts == 4])
print(paste("被全部4个模型选中的核心下调Hub代谢物数量:", length(down_core_hub)))

# 保存下调Hub代谢物
down_hub_df <- data.frame(
  Feature = names(down_feature_counts),
  Count = as.numeric(down_feature_counts),
  Regulation = "Down-regulated"
)
down_hub_df <- down_hub_df[order(down_hub_df$Count, decreasing = TRUE), ]
write.csv(down_hub_df, "Lipid_ML_Results/data/down_hub_metabolites.csv", row.names = FALSE)

# --------------------------------
# 8.3 合并所有Hub代谢物
# --------------------------------
print("合并上调和下调Hub代谢物...")

# 组合上调和下调hub代谢物
combined_hub_df <- rbind(up_hub_df, down_hub_df)
combined_hub_df <- combined_hub_df[order(combined_hub_df$Count, decreasing = TRUE), ]
write.csv(combined_hub_df, "Lipid_ML_Results/data/combined_hub_metabolites.csv", row.names = FALSE)

# 提取所有hub代谢物特征列表
all_hub_features <- c(up_hub_features, down_hub_features)
print(paste("总Hub代谢物数量:", length(all_hub_features)))

# 提取高置信度hub代谢物
high_confidence_hub_features <- c(up_high_confidence_hub, down_high_confidence_hub)
print(paste("高置信度Hub代谢物数量:", length(high_confidence_hub_features)))

# 提取核心hub代谢物
core_hub_features <- c(up_core_hub, down_core_hub)
print(paste("核心Hub代谢物数量:", length(core_hub_features)))

# ===============================
# 9. 创建可视化
# ===============================
print("创建可视化图表...")

# --------------------------------
# 9.1 创建韦恩图展示模型选择的特征交集
# --------------------------------
library(VennDiagram)
library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
library(ggpubr)

# 确保输出目录存在
dir.create("Lipid_ML_Results/figures", recursive = TRUE, showWarnings = FALSE)

# 选择4个模型的结果进行对比
# 对上调的脂质特征
up_selected_models <- list(
  "RF" = rf_up_selected,
  "Elastic Net" = elastic_up_selected,
  "SVM" = svm_up_selected,
  "XGBoost" = xgb_up_selected
)

# 对下调的脂质特征
down_selected_models <- list(
  "RF" = rf_down_selected,
  "Elastic Net" = elastic_down_selected,
  "SVM" = svm_down_selected,
  "XGBoost" = xgb_down_selected
)

# 设置颜色方案
venn_colors <- c(
  "RF" = "#E41A1C",          # 红色
  "Elastic Net" = "#377EB8", # 蓝色
  "SVM" = "#984EA3",         # 紫色
  "XGBoost" = "#4DAF4A"      # 绿色
)

# 创建上调脂质的韦恩图
venn.plot.up <- venn.diagram(
  x = up_selected_models,
  filename = NULL,  # 暂不保存，后面保存处理后的图
  category.names = names(up_selected_models),
  fill = venn_colors[names(up_selected_models)],
  alpha = 0.6,
  cat.col = venn_colors[names(up_selected_models)],
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.15,
  main = "Up-regulated Lipid Features Selected by Different Models",
  main.cex = 1.5,
  main.fontface = "bold",
  euler.d = TRUE,      # 使用欧拉图更精确显示关系
  scaled = TRUE,       # 按比例缩放圆的大小
  ext.text = FALSE,    # 不显示外部标签
  cat.pos = c(0, 0, 180, 180), # 调整标签位置
  cat.dist = c(0.1, 0.1, 0.1, 0.1) # 调整标签与圆的距离
)

# 创建下调脂质的韦恩图
venn.plot.down <- venn.diagram(
  x = down_selected_models,
  filename = NULL,  # 暂不保存，后面保存处理后的图
  category.names = names(down_selected_models),
  fill = venn_colors[names(down_selected_models)],
  alpha = 0.6,
  cat.col = venn_colors[names(down_selected_models)],
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.15,
  main = "Down-regulated Lipid Features Selected by Different Models",
  main.cex = 1.5,
  main.fontface = "bold",
  euler.d = TRUE,
  scaled = TRUE,
  ext.text = FALSE,
  cat.pos = c(0, 0, 180, 180),
  cat.dist = c(0.1, 0.1, 0.1, 0.1)
)

# 保存上调脂质的韦恩图
pdf("Lipid_ML_Results/figures/up_venn_diagram.pdf", width = 10, height = 8)
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# 保存下调脂质的韦恩图
pdf("Lipid_ML_Results/figures/down_venn_diagram.pdf", width = 10, height = 8)
grid.newpage()
grid.draw(venn.plot.down)
dev.off()


# --------------------------------
# 9.3 创建Hub代谢物热图
# --------------------------------
# 9.5 创建高置信度Hub代谢物表达热图(分组版)
# --------------------------------
# 为热图准备注释数据
annotation_col <- data.frame(
  HBP_Status = factor(diethabit_prolip$HBP, levels = c(0, 1), labels = c("Control", "Case"))
)
rownames(annotation_col) <- rownames(high_conf_hub_scaled)

# 创建分组热图 - 按HBP状态分组
# 找出control组样本和case组样本
control_samples <- which(annotation_col$HBP_Status == "Control")
case_samples <- which(annotation_col$HBP_Status == "Case")

# 对control组和case组的样本分别进行聚类
control_data <- high_conf_hub_scaled[control_samples, ]
case_data <- high_conf_hub_scaled[case_samples, ]

# 对control组进行层次聚类
control_dist <- dist(control_data)
control_hclust <- hclust(control_dist, method = "ward.D2")
control_order <- control_samples[control_hclust$order]

# 对case组进行层次聚类
case_dist <- dist(case_data)
case_hclust <- hclust(case_dist, method = "ward.D2")
case_order <- case_samples[case_hclust$order]

# 合并排序样本 - 先control后case，但各自按聚类顺序排列
ordered_idx <- c(control_order, case_order)

# 创建高置信度Hub代谢物表达矩阵
hub_matrix <- t(high_conf_hub_scaled[ordered_idx, ])

# 准备新的列注释
annotation_col_ordered <- data.frame(
  HBP_Status = annotation_col$HBP_Status[ordered_idx]
)
rownames(annotation_col_ordered) <- colnames(hub_matrix)

# 创建行注释 - 代谢物调控状态
regulation_status <- rep("Down-regulated", nrow(hub_matrix))
regulation_status[rownames(hub_matrix) %in% up_hub_features] <- "Up-regulated"

annotation_row <- data.frame(
  Regulation = factor(regulation_status)
)
rownames(annotation_row) <- rownames(hub_matrix)

# 创建按HBP分组的热图
hub_heatmap_by_hbp <- pheatmap(
  hub_matrix,
  cluster_rows = TRUE,      # 对代谢物(行)进行聚类
  cluster_cols = FALSE,     # 不对样本(列)进行聚类，我们已经手动聚类了
  scale = "row",            # 按行（代谢物）缩放
  annotation_col = annotation_col_ordered,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  main = "High Confidence Hub Lipid Metabolites by HBP Status",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  fontsize = 10,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-3, 3, length.out = 101),
  border_color = "grey80",
  gaps_col = length(control_order),  # 在对照组和病例组之间添加间隔
  clustering_method = "ward.D2",
  clustering_distance_rows = "correlation",
  annotation_names_col = TRUE,
  cellwidth = 5,
  cellheight = 12
)

# 保存按HBP分组的热图
pdf("Lipid_ML_Results/figures/hub_metabolites_by_hbp_groups.pdf", width = 12, height = 10)
print(hub_heatmap_by_hbp)

# 保存为PPT
# 使用eoffice保存ROC曲线图到PPT
topptx(hub_heatmap_by_hbp, 
       filename = "Lipid_ML_Results/figures/hub_metabolites_by_hbp_within_group_clustering.pptx",
       width = 7, 
       height = 6,
       append = FALSE,
       devsize = FALSE,
       units = "in")

# --------------------------------
# 9.3 创建箱线图展示Hub代谢物表达差异 (修订版)
# --------------------------------
# 从数据集中提取共有Hub代谢物数据
common_hub_data <- matrix(NA, nrow = nrow(diethabit_prolip), ncol = length(common_all_hub))
colnames(common_hub_data) <- common_all_hub
rownames(common_hub_data) <- rownames(diethabit_prolip)

# 填充数据
for (metabolite in common_all_hub) {
  if (metabolite %in% colnames(diethabit_prolip)) {
    common_hub_data[, metabolite] <- as.numeric(diethabit_prolip[, metabolite])
  }
}

# 准备箱线图数据
boxplot_data <- as.data.frame(common_hub_data)
boxplot_data$HBP_Status <- factor(diethabit_prolip$HBP, levels = c(0, 1), labels = c("control", "case"))

# 转换为长格式
boxplot_long <- reshape2::melt(boxplot_data, id.vars = "HBP_Status",
                               variable.name = "Metabolite", value.name = "Expression")

# 计算每个代谢物的p值
p_values <- list()
for(metabolite in unique(boxplot_long$Metabolite)) {
  subset_data <- boxplot_long[boxplot_long$Metabolite == metabolite, ]
  test_result <- t.test(Expression ~ HBP_Status, data = subset_data)
  p_values[[metabolite]] <- test_result$p.value
}

# 创建p值标记函数
get_significance <- function(p) {
  if(p < 0.0001) return("****")
  else if(p < 0.001) return("***")
  else if(p < 0.01) return("**")
  else if(p < 0.05) return("*")
  else return("ns")
}

# 为每个代谢物添加显著性标记
sig_data <- data.frame(
  Metabolite = names(p_values),
  p_value = unlist(p_values),
  significance = sapply(unlist(p_values), get_significance)
)

# 计算每个代谢物的数据范围，用于正确放置显著性标记
y_ranges <- boxplot_long %>%
  group_by(Metabolite) %>%
  summarize(
    min_val = min(Expression, na.rm = TRUE),
    max_val = max(Expression, na.rm = TRUE),
    range = max_val - min_val,
    # 设置标记位置（远离数据点但在图表范围内）
    mark_y = max_val + 0.1 * range
  )

# 合并显著性和Y轴位置信息
sig_data <- merge(sig_data, y_ranges, by = "Metabolite")

# 创建改进的箱线图 - 确保每个分面的Y轴根据数据自适应
p_improved <- ggplot(boxplot_long, aes(x = HBP_Status, y = Expression, fill = HBP_Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, aes(color = HBP_Status)) +
  # 不使用stat_compare_means，而是手动添加标记
  facet_wrap(~ Metabolite, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("control" = "#3498DB", "case" = "#E74C3C")) +
  scale_color_manual(values = c("control" = "#3498DB", "case" = "#E74C3C")) +
  labs(
    title = "Common Hub Metabolites Expression by HBP Status",
    x = "HBP Status",
    y = "Expression Level"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

# 手动添加显著性标记
# 遍历每个代谢物
for (i in 1:nrow(sig_data)) {
  metabolite <- sig_data$Metabolite[i]
  sig_mark <- sig_data$significance[i]
  y_pos <- sig_data$mark_y[i]
  
  # 只添加非NS标记
  if (sig_mark != "ns") {
    # 添加显著性标记和线段
    # 线段位置微调以获得良好视觉效果
    p_improved <- p_improved + 
      geom_text(
        data = data.frame(
          x = 1.5,
          y = y_pos,
          Metabolite = metabolite,
          label = sig_mark
        ),
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 3.5
      ) +
      geom_segment(
        data = data.frame(
          x = 1,
          xend = 2,
          y = y_pos * 0.95,
          yend = y_pos * 0.95,
          Metabolite = metabolite
        ),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        size = 0.5
      )
  }
}

# 保存图表
ggsave("Lipid_ML_Results/figures/common_hub_boxplots_improved.pdf", p_improved, width = 12, height = 10)
# 使用eoffice保存到PPT
topptx(p_improved, 
       filename = "Lipid_ML_Results/figures/p_improved.pptx",
       width = 8, 
       height = 8,
       append = FALSE,
       devsize = FALSE,
       units = "in")
# --------------------------------
# 9.3 创建不分面的核心Hub代谢物表达差异箱线图
# --------------------------------
library(ggplot2)
library(dplyr)
library(reshape2)

# 定义核心Hub代谢物
common_up_hub <- c("PE(O/16:0/18:1)")
common_down_hub <- c("TAG(56:6)_FA22:5", "TAG(50:2)_FA14:0", "TAG(52:2)_FA14:0", "TAG(58:6)_FA18:0")
common_all_hub <- c(common_up_hub, common_down_hub)

# 从数据集中提取核心Hub代谢物数据
hub_data <- as.data.frame(hub_lipids_scaled)
hub_data$HBP_Status <- factor(diethabit_prolip$HBP, levels = c(0, 1), labels = c("Control", "Case"))

# 只保留核心Hub代谢物
core_columns <- c("HBP_Status", common_all_hub)
core_data_filtered <- hub_data[, core_columns]

# 转换为长格式用于ggplot2
core_data_long <- reshape2::melt(core_data_filtered, id.vars = "HBP_Status",
                                 variable.name = "Metabolite", value.name = "Expression")

# 添加调控信息（用于颜色区分，不用于分面）
core_data_long$Regulation <- "Down-regulated"
core_data_long$Regulation[core_data_long$Metabolite %in% common_up_hub] <- "Up-regulated"

# 将Regulation转换为因子
core_data_long$Regulation <- factor(core_data_long$Regulation, 
                                    levels = c("Up-regulated", "Down-regulated"))

# 计算每个代谢物的p值
p_values <- list()
for(metabolite in unique(core_data_long$Metabolite)) {
  subset_data <- core_data_long[core_data_long$Metabolite == metabolite, ]
  test_result <- t.test(Expression ~ HBP_Status, data = subset_data)
  p_values[[as.character(metabolite)]] <- test_result$p.value
}

# 创建p值标记函数
get_significance <- function(p) {
  if(p < 0.0001) return("****")
  else if(p < 0.001) return("***")
  else if(p < 0.01) return("**")
  else if(p < 0.05) return("*")
  else return("ns")
}

# 为每个代谢物添加显著性标记
sig_data <- data.frame(
  Metabolite = names(p_values),
  p_value = unlist(p_values),
  significance = sapply(unlist(p_values), get_significance),
  stringsAsFactors = FALSE
)

# 创建不分面的箱线图
# 可以选择性地按调控方向重新排序代谢物
# 例如，先显示所有上调代谢物，再显示所有下调代谢物
core_data_long$Metabolite <- factor(core_data_long$Metabolite, 
                                    levels = common_all_hub)

# 创建单个箱线图
core_boxplot <- ggplot(core_data_long, aes(x = Metabolite, y = Expression, fill = HBP_Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = HBP_Status), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
              alpha = 0.6, size = 0.8) +
  scale_fill_manual(values = c("Control" = "#3498DB", "Case" = "#E74C3C"), name = "HBP Status") +
  scale_color_manual(values = c("Control" = "#3498DB", "Case" = "#E74C3C"), name = "HBP Status") +
  # 使用形状或边框颜色来区分上调和下调代谢物
  geom_point(data = core_data_long[core_data_long$Regulation == "Up-regulated",], 
             aes(x = Metabolite, y = -2.5), shape = 24, fill = "darkgreen", size = 3, alpha = 0.7) +
  labs(
    title = "Core Hub Lipid Metabolites",
    subtitle = "Selected by all 4 machine learning models (n=5)",
    x = "",
    y = "Standardized Expression"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90")
  )

# 添加显著性标记
for(i in 1:nrow(sig_data)) {
  metabolite <- sig_data$Metabolite[i]
  sig_mark <- sig_data$significance[i]
  
  # 只为显著的代谢物添加标记
  if(sig_mark != "ns") {
    # 计算该代谢物的y轴最大值
    max_y <- max(core_data_long$Expression[core_data_long$Metabolite == metabolite], na.rm = TRUE)
    
    # 添加显著性标记
    core_boxplot <- core_boxplot +
      annotate("text", 
               x = metabolite, 
               y = max_y + 0.6, 
               label = sig_mark, 
               size = 3.5)
  }
}

# 添加上调和下调的标注
# 创建一个简单的图例说明
core_boxplot <- core_boxplot +
  annotate("text", x = 0.8, y = -2.5, label = "Up-regulated", color = "darkgreen", size = 3.5, hjust = 1) +
  # 添加图例边框
  annotate("segment", x = 0.5, xend = 1.1, y = -2.7, yend = -2.7, color = "black", size = 0.5) +
  annotate("segment", x = 0.5, xend = 1.1, y = -2.3, yend = -2.3, color = "black", size = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = -2.7, yend = -2.3, color = "black", size = 0.5) +
  annotate("segment", x = 1.1, xend = 1.1, y = -2.7, yend = -2.3, color = "black", size = 0.5)

# 保存核心Hub代谢物箱线图
ggsave("Lipid_ML_Results/figures/core_hub_boxplot_single.pdf", core_boxplot, width = 12, height = 8)
ggsave("Lipid_ML_Results/figures/core_hub_boxplot_single.png", core_boxplot, width = 12, height = 8, dpi = 300)

# --------------------------------
# 同样为高置信度Hub代谢物创建单个箱线图
# --------------------------------
# 定义高置信度Hub代谢物
high_confidence_up_hub <- up_high_confidence_hub  # 5个上调高置信度代谢物
high_confidence_down_hub <- down_high_confidence_hub  # 12个下调高置信度代谢物
high_confidence_all_hub <- c(high_confidence_up_hub, high_confidence_down_hub)

# 只保留高置信度Hub代谢物
high_conf_columns <- c("HBP_Status", high_confidence_all_hub)
high_conf_data_filtered <- hub_data[, high_conf_columns]

# 转换为长格式用于ggplot2
high_conf_data_long <- reshape2::melt(high_conf_data_filtered, id.vars = "HBP_Status",
                                      variable.name = "Metabolite", value.name = "Expression")

# 添加调控信息
high_conf_data_long$Regulation <- "Down-regulated"
high_conf_data_long$Regulation[high_conf_data_long$Metabolite %in% high_confidence_up_hub] <- "Up-regulated"

# 将Regulation转换为因子
high_conf_data_long$Regulation <- factor(high_conf_data_long$Regulation, 
                                         levels = c("Up-regulated", "Down-regulated"))

# 为高置信度Hub代谢物计算p值
high_conf_p_values <- list()
for(metabolite in unique(high_conf_data_long$Metabolite)) {
  subset_data <- high_conf_data_long[high_conf_data_long$Metabolite == metabolite, ]
  test_result <- t.test(Expression ~ HBP_Status, data = subset_data)
  high_conf_p_values[[as.character(metabolite)]] <- test_result$p.value
}

# 为高置信度Hub代谢物添加显著性标记
high_conf_sig_data <- data.frame(
  Metabolite = names(high_conf_p_values),
  p_value = unlist(high_conf_p_values),
  significance = sapply(unlist(high_conf_p_values), get_significance),
  stringsAsFactors = FALSE
)

# 创建高置信度代谢物的单个箱线图
high_conf_boxplot <- ggplot(high_conf_data_long, aes(x = Metabolite, y = Expression, fill = HBP_Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = HBP_Status), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
              alpha = 0.6, size = 0.8) +
  scale_fill_manual(values = c("Control" = "#3498DB", "Case" = "#E74C3C"), name = "HBP Status") +
  scale_color_manual(values = c("Control" = "#3498DB", "Case" = "#E74C3C"), name = "HBP Status") +
  # 使用形状标记上调代谢物
  geom_point(data = high_conf_data_long[high_conf_data_long$Regulation == "Up-regulated",], 
             aes(x = Metabolite, y = -2.5), shape = 24, fill = "darkgreen", size = 3, alpha = 0.7) +
  labs(
    title = "High Confidence Hub Lipid Metabolites",
    subtitle = "Selected by at least 3 different machine learning models (n=17)",
    x = "",
    y = "Standardized Expression"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90")
  )

# 添加高置信度Hub代谢物的显著性标记
for(i in 1:nrow(high_conf_sig_data)) {
  metabolite <- high_conf_sig_data$Metabolite[i]
  sig_mark <- high_conf_sig_data$significance[i]
  
  # 只为显著的代谢物添加标记
  if(sig_mark != "ns") {
    # 计算该代谢物的y轴最大值
    max_y <- max(high_conf_data_long$Expression[high_conf_data_long$Metabolite == metabolite], na.rm = TRUE)
    
    # 添加显著性标记
    high_conf_boxplot <- high_conf_boxplot +
      annotate("text", 
               x = metabolite, 
               y = max_y + 0.6, 
               label = sig_mark, 
               size = 3.5)
  }
}

# 添加上调和下调的标注
high_conf_boxplot <- high_conf_boxplot +
  annotate("text", x = 0.8, y = -2.5, label = "Up-regulated", color = "darkgreen", size = 3.5, hjust = 1) +
  # 添加图例边框
  annotate("segment", x = 0.5, xend = 1.1, y = -2.7, yend = -2.7, color = "black", size = 0.5) +
  annotate("segment", x = 0.5, xend = 1.1, y = -2.3, yend = -2.3, color = "black", size = 0.5) +
  annotate("segment", x = 0.5, xend = 0.5, y = -2.7, yend = -2.3, color = "black", size = 0.5) +
  annotate("segment", x = 1.1, xend = 1.1, y = -2.7, yend = -2.3, color = "black", size = 0.5)

# 保存高置信度Hub代谢物箱线图
ggsave("Lipid_ML_Results/figures/high_confidence_hub_boxplot_single.pdf", high_conf_boxplot, width = 16, height = 9)
ggsave("Lipid_ML_Results/figures/high_confidence_hub_boxplot_single.png", high_conf_boxplot, width = 16, height = 9, dpi = 300)
# --------------------------------
# 9.4 创建顶部特征的棒棒糖图
# --------------------------------
# 准备各模型的重要性数据
# 随机森林
rf_up_top10 <- head(rf_up_importance_df, 15)
rf_up_top10$Model <- "Random Forest"
rf_up_top10$Regulation <- "Up-regulated"

rf_down_top10 <- head(rf_down_importance_df, 15)
rf_down_top10$Model <- "Random Forest"
rf_down_top10$Regulation <- "Down-regulated"

# 弹性网络 - 使用AbsCoef作为重要性指标
elastic_up_top10 <- head(elastic_up_importance_df, 15)
elastic_up_top10$Model <- "Elastic Net"
elastic_up_top10$Regulation <- "Up-regulated"
elastic_up_top10$Importance <- elastic_up_top10$AbsCoef

elastic_down_top10 <- head(elastic_down_importance_df, 15)
elastic_down_top10$Model <- "Elastic Net"
elastic_down_top10$Regulation <- "Down-regulated"
elastic_down_top10$Importance <- elastic_down_top10$AbsCoef

# XGBoost
xgb_up_top10 <- head(xgb_up_importance_df, 15)
xgb_up_top10$Model <- "XGBoost"
xgb_up_top10$Regulation <- "Up-regulated"

xgb_down_top10 <- head(xgb_down_importance_df, 15)
xgb_down_top10$Model <- "XGBoost"
xgb_down_top10$Regulation <- "Down-regulated"

# SVM-RFE - 使用排名作为重要性
svm_up_top10 <- head(svm_up_importance, 15)
svm_up_top10$Model <- "SVM-RFE"
svm_up_top10$Regulation <- "Up-regulated"

svm_down_top10 <- head(svm_down_importance, 15)
svm_down_top10$Model <- "SVM-RFE"
svm_down_top10$Regulation <- "Down-regulated"

# 合并所有数据
all_models_top10 <- rbind(
  rf_up_top10[, c("Feature", "Importance", "Model", "Regulation")],
  rf_down_top10[, c("Feature", "Importance", "Model", "Regulation")],
  elastic_up_top10[, c("Feature", "Importance", "Model", "Regulation")],
  elastic_down_top10[, c("Feature", "Importance", "Model", "Regulation")],
  xgb_up_top10[, c("Feature", "Importance", "Model", "Regulation")],
  xgb_down_top10[, c("Feature", "Importance", "Model", "Regulation")],
  svm_up_top10[, c("Feature", "Importance", "Model", "Regulation")],
  svm_down_top10[, c("Feature", "Importance", "Model", "Regulation")]
)

# 为了使图形美观，对每个模型的重要性分数进行标准化
all_models_top10 <- all_models_top10 %>%
  group_by(Model, Regulation) %>%
  mutate(
    Importance_Scaled = (Importance - min(Importance)) / (max(Importance) - min(Importance)) * 100
  ) %>%
  ungroup()

# 设置模型顺序
all_models_top10$Model <- factor(all_models_top10$Model, 
                                 levels = c("Random Forest", "Elastic Net", "SVM-RFE", "XGBoost"))

# 设置颜色方案
model_colors <- c(
  "Random Forest" = "#E41A1C", 
  "Elastic Net" = "#377EB8", 
  "SVM-RFE" = "#984EA3", 
  "XGBoost" = "#4DAF4A"
)

# 创建顶部特征棒棒糖图
p_lollipop <- ggplot(all_models_top10, 
                     aes(x = reorder(Feature, Importance_Scaled), 
                         y = Importance_Scaled, 
                         color = Model)) +
  geom_segment(aes(x = reorder(Feature, Importance_Scaled), 
                   xend = reorder(Feature, Importance_Scaled), 
                   y = 0, 
                   yend = Importance_Scaled), 
               color = "gray50", size = 0.7) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = model_colors) +
  coord_flip() +
  facet_grid(Regulation ~ Model, scales = "free_y") +
  labs(
    title = "Top 10 Lipid Features Selected by Different Machine Learning Models",
    x = "Metabolite Feature", 
    y = "Relative Importance (%)", 
    color = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "#F0F0F0"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

# 保存棒棒糖图
ggsave("Lipid_ML_Results/figures/models_feature_selection_lollipop.pdf", p_lollipop, width = 16, height = 12)

