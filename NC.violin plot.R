
# 清空环境
rm(list = ls())
gc()
# 加载必要的包
required_packages <- c(
  "data.table", "dplyr", "ggplot2", "rstatix", "ggpubr", "ggh4x",
  "officer", "rvg", "gridExtra", "scales", "RColorBrewer"
)

# 检查并安装缺失的包
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# 加载包
lapply(required_packages, library, character.only = TRUE)
loaded_packages <- search()


# =============================================================================
# 1. 数据读取和预处理（使用您现有的代码）TREData
# =============================================================================

cat("\n=== 读取和预处理数据 ===\n")

# 读取原始数据
data <- read_excel("TREData.xlsx")

# 数据预处理
data <- data %>%
  mutate(
    # 主要饮食时间分类变量
    eattime01 = factor(eattime01, 
                       levels = c(0, 1), 
                       labels = c("≤12.5h", ">12.5h")),
    
    eatt301 = factor(eatt301, 
                     levels = c(0, 1), 
                     labels = c("≤20:00", ">20:00")),
    
    # 协变量
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
    urban = factor(urban, levels = c(1, 2), labels = c("Rural", "Urban")),
    age12 = factor(age12, levels = c(1, 2), labels = c("<12 years", "≥12 years")),
    fhistory = as.factor(fhistory),
    MARRIAGE = as.factor(MARRIAGE),
    FEDU = as.factor(FEDU),
    MEDU = as.factor(MEDU),
    FOCUP = as.factor(FOCUP),
    MOCUP = as.factor(MOCUP),
    INCOME = as.factor(INCOME),
    
    # 确保心血管指标为数值型
    across(c(SBP_mean, DBP_mean, R_pwv, L_pwv, cIMT_L, cIMT_R, 
             LVDd, LVDs, EF, FS, CHOL, TG, HDL, LDL, GLU), as.numeric),
    
    # 确保数值型协变量为数值型
    across(c(height_zscore, WHtR, Energy, chei_score_raw, 
             protein_energy_ratio, 
             total_met), as.numeric)
  )

# 定义变量
all_cv_vars <- c("SBP_mean", "DBP_mean", "R_pwv", "L_pwv", 
                 "cIMT_L", "cIMT_R", "EF", "FS", "LVDd", "LVDs",
                 "CHOL", "TG", "HDL", "LDL", "GLU")

dietary_vars <- c("eattime01", "eatt301")

covariates <- c("age12", "height_zscore", "sex", "WHtR", "urban", "fhistory", 
                "MARRIAGE", "FEDU", "MEDU", "FOCUP", "MOCUP", "INCOME", 
                "Energy", "chei_score_raw", "protein_energy_ratio", "total_met")

# =============================================================================
# 2. 协方差分析函数（使用您现有的函数）
# =============================================================================

perform_ancova <- function(data, outcome_var, group_var, covariates) {
  
  tryCatch({
    # 只排除自变量和因变量缺失的样本
    essential_vars <- c(outcome_var, group_var)
    all_vars <- c(essential_vars, covariates[covariates %in% names(data)])
    
    analysis_data <- data %>%
      select(all_of(all_vars)) %>%
      filter(!is.na(!!sym(outcome_var)) & !is.na(!!sym(group_var)))
    
    if(nrow(analysis_data) < 20) {
      return(list(
        success = FALSE,
        message = paste("样本量不足:", nrow(analysis_data)),
        overall_p = NA,
        effect_size = NA,
        group_means = NA,
        group_se = NA
      ))
    }
    
    cat("   - 分析样本量:", nrow(analysis_data), "(原始:", nrow(data), ")\n")
    
    # 检查协变量可用性和变异性
    available_covs <- c()
    for(cov in covariates) {
      if(cov %in% names(analysis_data)) {
        non_missing_values <- analysis_data[[cov]][!is.na(analysis_data[[cov]])]
        
        if(length(non_missing_values) >= 10) {
          if(is.numeric(analysis_data[[cov]])) {
            if(var(non_missing_values, na.rm = TRUE) > 0) {
              available_covs <- c(available_covs, cov)
            }
          } else if(is.factor(analysis_data[[cov]])) {
            if(length(unique(non_missing_values)) > 1) {
              available_covs <- c(available_covs, cov)
            }
          }
        }
      }
    }
    
    # 构建ANCOVA模型公式
    if(length(available_covs) > 0) {
      formula_str <- paste(outcome_var, "~", group_var, "+", paste(available_covs, collapse = " + "))
      cat("   - 使用协变量:", length(available_covs), "个\n")
    } else {
      formula_str <- paste(outcome_var, "~", group_var)
      cat("   - 无可用协变量，执行简单方差分析\n")
    }
    
    # 拟合模型
    model <- aov(as.formula(formula_str), data = analysis_data)
    actual_n <- length(model$residuals)
    
    # 总体F检验
    model_summary <- summary(model)
    group_p_value <- model_summary[[1]][group_var, "Pr(>F)"]
    
    # 效应量计算
    if("rstatix" %in% loaded_packages) {
      eta_squared <- tryCatch({
        result <- rstatix::eta_squared(model)
        if(nrow(result) > 0) result$ges[1] else NA
      }, error = function(e) NA)
    } else {
      eta_squared <- NA
    }
    
    # 计算组均值和标准误
    group_stats <- analysis_data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        mean_val = mean(!!sym(outcome_var), na.rm = TRUE),
        se_val = sd(!!sym(outcome_var), na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = 'drop'
      )
    
    return(list(
      success = TRUE,
      model = model,
      overall_p = group_p_value,
      effect_size = eta_squared,
      n = nrow(analysis_data),
      actual_n = actual_n,
      formula = formula_str,
      group_means = setNames(group_stats$mean_val, group_stats[[group_var]]),
      group_se = setNames(group_stats$se_val, group_stats[[group_var]]),
      available_covariates = available_covs,
      covariate_missing_rate = sapply(available_covs, function(x) {
        round(sum(is.na(analysis_data[[x]])) / nrow(analysis_data) * 100, 1)
      }),
      analysis_data = analysis_data
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      message = paste("ANCOVA失败:", e$message),
      overall_p = NA,
      effect_size = NA
    ))
  })
}

# =============================================================================
# 3. 美化变量名函数
# =============================================================================

beautify_variable_name <- function(y_var) {
  case_when(
    y_var == "SBP_mean" ~ "SBP (mmHg)",
    y_var == "DBP_mean" ~ "DBP (mmHg)",
    y_var == "R_pwv" ~ "Right baPWV (m/s)",
    y_var == "L_pwv" ~ "Left baPWV (m/s)",
    y_var == "cIMT_L" ~ "Left cIMT (mm)",
    y_var == "cIMT_R" ~ "Right cIMT (mm)",
    y_var == "EF" ~ "EF (%)",
    y_var == "FS" ~ "FS (%)",
    y_var == "LVDd" ~ "LVDd (mm)",
    y_var == "LVDs" ~ "LVDs (mm)",
    y_var == "CHOL" ~ "TC (mg/dL)",
    y_var == "TG" ~ "TG (mg/dL)",
    y_var == "HDL" ~ "HDL-C (mg/dL)",
    y_var == "LDL" ~ "LDL-C (mg/dL)",
    y_var == "GLU" ~ "GLU (mg/dL)",
    TRUE ~ y_var
  )
}

beautify_group_name <- function(group_var) {
  case_when(
    group_var == "eattime01" ~ "Eating Duration",
    group_var == "eatt301" ~ "Late Eating Time",
    TRUE ~ group_var
  )
}

# =============================================================================
# 4. 可视化函数
# =============================================================================

# 创建单个子图的函数（适用于拼接）
create_cardiovascular_subplot <- function(data, outcome_var, group_var, ancova_result = NULL) {
  
  # 美化变量名
  y_label <- beautify_variable_name(outcome_var)
  
  # 为每个饮食时间变量设置不同的颜色
  if(group_var == "eattime01") {
    fill_colors <- c("#3B9AB2", "#E6A0C4")  # 蓝色和粉色
  } else {
    fill_colors <- c("#7294D4", "#FCAE12")  # 淡蓝色和橙色
  }
  
  # 准备数据进行分析
  plot_data <- data %>%
    filter(!is.na(!!sym(outcome_var)) & !is.na(!!sym(group_var))) %>%
    select(!!sym(outcome_var), !!sym(group_var))
  
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # 进行统计检验
  stat_test <- NULL
  tryCatch({
    stat_test <- plot_data %>%
      wilcox_test(as.formula(paste(outcome_var, "~", group_var))) %>%
      add_significance() %>%
      add_xy_position(x = group_var, dodge = 0.8)
  }, error = function(e) {
    cat("统计检验失败:", e$message, "\n")
  })
  
  # 获取ANCOVA的p值和显著性
  ancova_sig <- "ns"
  if(!is.null(ancova_result) && ancova_result$success) {
    ancova_p <- ancova_result$overall_p
    ancova_sig <- case_when(
      ancova_p < 0.001 ~ "***",
      ancova_p < 0.01 ~ "**", 
      ancova_p < 0.05 ~ "*",
      ancova_p < 0.1 ~ "†",
      TRUE ~ "ns"
    )
  }
  
  # 创建基础图形 - 调整为适合拼接的子图
  p <- plot_data %>%
    ggplot(aes(x = !!sym(group_var), y = !!sym(outcome_var))) +
    
    # 分布图表
    geom_violin(aes(fill = !!sym(group_var)), 
                trim = FALSE, 
                show.legend = FALSE,
                alpha = 0.7) +
    geom_boxplot(width = 0.3, 
                 outlier.shape = NA, 
                 staplewidth = 0.5,
                 alpha = 0.8,
                 size = 0.8) +
    
    # 统计元素
    stat_summary(fun = mean, 
                 geom = "point",
                 color = "black", 
                 fill = "#F98400",
                 shape = 23, 
                 size = 4, 
                 show.legend = FALSE) +
    
    # 填充颜色
    scale_fill_manual(values = fill_colors) +
    
    # 标签设置 - 只显示Y轴标签，标题包含显著性信息
    labs(
      x = "",  # 移除X轴标签以节省空间
      y = y_label,
      title = paste0(y_label, " (", ancova_sig, ")")
    ) +
    
    # 主题设置 - 适合拼接的大字号设置
    theme_bw() +
    theme(
      # 文本设置 - 增大字号
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                                 color = "black", size = 16, face = "bold"),
      axis.text.y = element_text(color = "black", size = 14, face = "bold"),
      axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 8)),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      
      # 背景和网格
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "grey80"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 边距 - 减小边距以适应拼接
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
      
      # 坐标轴线
      axis.line = element_line(color = "black", size = 0.8),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2)
    )
  
  # 添加统计显著性标记
  if(!is.null(stat_test)) {
    tryCatch({
      p <- p + stat_pvalue_manual(stat_test,
                                  label = "p.signif",
                                  hide.ns = TRUE,
                                  tip.length = 0.01,
                                  label.size = 6,
                                  color = "black",
                                  y.position = max(plot_data[[outcome_var]], na.rm = TRUE) * 1.08)
    }, error = function(e) {
      cat("添加显著性标记失败:", e$message, "\n")
    })
  }
  
  return(p)
}

# 创建拼接大图的函数
create_combined_cardiovascular_plot <- function(data, group_var, all_cv_vars, all_ancova_results) {
  
  cat("创建", beautify_group_name(group_var), "的拼接大图...\n")
  
  # 创建所有子图
  subplot_list <- list()
  
  for(cv_var in all_cv_vars) {
    result_key <- paste(cv_var, group_var, sep = "_vs_")
    ancova_result <- all_ancova_results[[result_key]]
    
    if(!is.null(ancova_result) && ancova_result$success) {
      subplot <- create_cardiovascular_subplot(data, cv_var, group_var, ancova_result)
      if(!is.null(subplot)) {
        subplot_list[[cv_var]] <- subplot
        cat("✓ 子图创建成功:", beautify_variable_name(cv_var), "\n")
      }
    }
  }
  
  if(length(subplot_list) == 0) {
    cat("✗ 没有可用的子图\n")
    return(NULL)
  }
  
  # 计算网格布局 - 3列5行或4列4行，根据子图数量调整
  n_plots <- length(subplot_list)
  if(n_plots <= 12) {
    ncol <- 3
    nrow <- ceiling(n_plots / ncol)
  } else {
    ncol <- 4
    nrow <- ceiling(n_plots / ncol)
  }
  
  # 创建主标题
  main_title <- paste("Cardiovascular Risk Factors vs", beautify_group_name(group_var))
  subtitle <- paste("ANCOVA Analysis Results (n =", nrow(data), ")")
  
  # 使用cowplot创建拼接图
  # 首先创建标题
  title_plot <- ggplot() + 
    labs(title = main_title, subtitle = subtitle) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 28, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 22, color = "gray40", margin = margin(b = 20))
    )
  
  # 创建网格拼接图
  combined_plot <- plot_grid(plotlist = subplot_list, 
                             ncol = ncol, 
                             nrow = nrow,
                             align = "hv",
                             axis = "tblr")
  
  # 添加组标签
  group_label <- beautify_group_name(group_var)
  if(group_var == "eattime01") {
    x_axis_labels <- c("≤12.5h", ">12.5h")
  } else {
    x_axis_labels <- c("≤20:00", ">20:00")
  }
  
  # 创建X轴标签
  x_label_plot <- ggplot() + 
    labs(title = paste(group_label, ":", paste(x_axis_labels, collapse = " vs "))) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(t = 15))
    )
  
  # 最终拼接：标题 + 图表 + X轴标签
  final_plot <- plot_grid(
    title_plot,
    combined_plot,
    x_label_plot,
    ncol = 1,
    rel_heights = c(0.12, 0.8, 0.08)
  )
  
  cat("✓ 拼接大图创建成功，包含", length(subplot_list), "个子图\n")
  
  return(final_plot)
}

# =============================================================================
# 5. PPT输出函数
# =============================================================================

# PPT输出函数 - 适应大图
create_enhanced_ppt_output <- function(plot_obj, filename, width = 20, height = 16, 
                                       units = "in", plot_type = "CardiovascularPlot") {
  
  tryCatch({
    # 创建PowerPoint对象
    ppt <- read_pptx()
    
    # 添加标题页
    ppt <- ppt %>%
      add_slide(layout = "Title Slide", master = "Office Theme") %>%
      ph_with(value = "Cardiovascular Risk Factors Analysis", 
              location = ph_location_type(type = "ctrTitle")) %>%
      ph_with(value = paste("Generated on", Sys.Date()), 
              location = ph_location_type(type = "subTitle"))
    
    # 添加图表页 - 使用空白布局以最大化图表空间
    ppt <- ppt %>%
      add_slide(layout = "Blank", master = "Office Theme") %>%
      ph_with(value = dml(ggobj = plot_obj, width = width*0.8, height = height*0.8), 
              location = ph_location(left = 0.5, top = 0.5, width = width*0.8, height = height*0.8))
    
    # 保存PPT文件
    ppt_filename <- paste0(filename, ".pptx")
    print(ppt, target = ppt_filename)
    
    # 保存高分辨率图片文件
    png_filename <- paste0(filename, ".png")
    ggsave(filename = png_filename, 
           plot = plot_obj, 
           width = width, 
           height = height, 
           dpi = 300, 
           device = 'png', 
           bg = '#FFFFFF')
    
    # 同时保存PDF版本（用于高质量打印）
    pdf_filename <- paste0(filename, ".pdf")
    ggsave(filename = pdf_filename, 
           plot = plot_obj, 
           width = width, 
           height = height, 
           device = 'pdf', 
           bg = '#FFFFFF')
    
    cat("✓ 已创建文件:", basename(ppt_filename), ",", basename(png_filename), "和", basename(pdf_filename), "\n")
    
    return(list(
      ppt_file = ppt_filename,
      png_file = png_filename,
      pdf_file = pdf_filename
    ))
    
  }, error = function(e) {
    cat("PPT创建失败:", e$message, "\n")
    return(NULL)
  })
}

# =============================================================================
# 6. 主分析流程 - 修改为创建拼接大图
# =============================================================================

cat("\n=== 开始心血管指标分析 ===\n")

# 创建输出目录
output_dir <- "cardiovascular_analysis_results"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 存储结果
all_ancova_results <- list()
created_files <- c()

# =============================================================================
# 第一步：执行所有协方差分析
# =============================================================================

cat("\n=== 第一步：执行所有协方差分析 ===\n")

for(cv_var in all_cv_vars) {
  for(diet_var in dietary_vars) {
    
    cat(paste("分析:", beautify_variable_name(cv_var), "vs", beautify_group_name(diet_var)), "\n")
    
    # 执行协方差分析
    ancova_result <- perform_ancova(data, cv_var, diet_var, covariates)
    
    if(ancova_result$success) {
      cat("✓ ANCOVA分析成功 - P值:", sprintf("%.4f", ancova_result$overall_p), "\n")
    } else {
      cat("✗ ANCOVA分析失败:", ancova_result$message, "\n")
    }
    
    # 存储结果
    result_key <- paste(cv_var, diet_var, sep = "_vs_")
    all_ancova_results[[result_key]] <- ancova_result
  }
}

# =============================================================================
# 第二步：创建拼接大图
# =============================================================================

cat("\n=== 第二步：创建拼接大图 ===\n")

# 为每个饮食时间变量创建拼接大图
for(diet_var in dietary_vars) {
  
  
  cat(paste("创建", beautify_group_name(diet_var), "的拼接大图"), "\n")
  
  
  # 创建拼接大图
  combined_plot <- create_combined_cardiovascular_plot(data, diet_var, all_cv_vars, all_ancova_results)
  
  if(is.null(combined_plot)) {
    cat("✗ 拼接大图创建失败\n")
    next
  }
  
  # 输出PPT和图像文件
  cat("输出文件...\n")
  filename_base <- file.path(output_dir, paste0("Combined_CV_Analysis_", diet_var))
  
  ppt_result <- create_enhanced_ppt_output(
    plot_obj = combined_plot,
    filename = filename_base,
    width = 24,  # 增大宽度以适应多个子图
    height = 18, # 增大高度以适应多个子图
    units = "in",
    plot_type = paste("Combined Cardiovascular Analysis -", beautify_group_name(diet_var))
  )
  
  if(!is.null(ppt_result)) {
    created_files <- c(created_files, 
                       basename(ppt_result$ppt_file), 
                       basename(ppt_result$png_file),
                       basename(ppt_result$pdf_file))
  }
  
  cat("✓", beautify_group_name(diet_var), "拼接大图创建完成\n")
}

# =============================================================================
# 7. 创建汇总报告
# =============================================================================

cat("\n=== 创建汇总报告 ===\n")

# 创建结果汇总表
summary_results <- data.frame(
  Cardiovascular_Variable = character(),
  Dietary_Variable = character(), 
  Sample_Size = numeric(),
  P_Value = numeric(),
  Effect_Size = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

for(result_name in names(all_results)) {
  result <- all_results[[result_name]]
  if(result$ancova_result$success) {
    significance <- case_when(
      result$ancova_result$overall_p < 0.001 ~ "***",
      result$ancova_result$overall_p < 0.01 ~ "**",
      result$ancova_result$overall_p < 0.05 ~ "*",
      result$ancova_result$overall_p < 0.1 ~ "†",
      TRUE ~ "ns"
    )
    
    summary_results <- rbind(summary_results, data.frame(
      Cardiovascular_Variable = beautify_variable_name(result$cardiovascular_var),
      Dietary_Variable = beautify_group_name(result$dietary_var),
      Sample_Size = result$ancova_result$n,
      P_Value = result$ancova_result$overall_p,
      Effect_Size = result$ancova_result$effect_size,
      Significance = significance,
      stringsAsFactors = FALSE
    ))
  }
}

# 保存汇总结果
write.csv(summary_results, file.path(output_dir, "summary_results1.csv"), row.names = FALSE)
created_files <- c(created_files, "summary_results1.csv")

# 输出最终汇总信息

cat("分析完成汇总\n")

cat("✓ 总共分析了", nrow(summary_results), "个变量组合\n")
cat("✓ 显著结果 (p<0.05):", sum(summary_results$P_Value < 0.05, na.rm = TRUE), "个\n")
cat("✓ 创建的文件:", length(created_files), "个\n")
cat("✓ 输出目录:", output_dir, "\n")

# 显示最显著的结果
if(nrow(summary_results) > 0) {
  significant_results <- summary_results[summary_results$P_Value < 0.05, ]
  if(nrow(significant_results) > 0) {
    cat("\n最显著的关联 (p<0.05):\n")
    print(significant_results[order(significant_results$P_Value), ])
  }
}

cat("\n=== 分析完成 ===\n")




# =============================================================================
# 7. 创建汇总报告 (修改版本 - 包含均值±SE)
# =============================================================================

cat("\n=== 创建汇总报告 ===\n")

# 创建结果汇总表 - 增加均值和标准误列
summary_results <- data.frame(
  Cardiovascular_Variable = character(),
  Dietary_Variable = character(), 
  Sample_Size = numeric(),
  Group1_Name = character(),
  Group1_Mean = numeric(),
  Group1_SE = numeric(),
  Group1_Mean_SE = character(),
  Group2_Name = character(),
  Group2_Mean = numeric(),
  Group2_SE = numeric(),
  Group2_Mean_SE = character(),
  P_Value = numeric(),
  Effect_Size = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

# 修正变量名错误：应该使用all_ancova_results而不是all_results
for(cv_var in all_cv_vars) {
  for(diet_var in dietary_vars) {
    
    result_key <- paste(cv_var, diet_var, sep = "_vs_")
    ancova_result <- all_ancova_results[[result_key]]
    
    if(!is.null(ancova_result) && ancova_result$success) {
      
      # 计算显著性
      significance <- case_when(
        ancova_result$overall_p < 0.001 ~ "***",
        ancova_result$overall_p < 0.01 ~ "**",
        ancova_result$overall_p < 0.05 ~ "*",
        ancova_result$overall_p < 0.1 ~ "†",
        TRUE ~ "ns"
      )
      
      # 获取组名和对应的均值、标准误
      group_names <- names(ancova_result$group_means)
      group_means <- ancova_result$group_means
      group_se <- ancova_result$group_se
      
      # 确保有两个组的数据
      if(length(group_names) >= 2) {
        
        # 格式化均值±SE
        group1_mean_se <- sprintf("%.2f ± %.2f", group_means[1], group_se[1])
        group2_mean_se <- sprintf("%.2f ± %.2f", group_means[2], group_se[2])
        
        # 添加到汇总结果
        summary_results <- rbind(summary_results, data.frame(
          Cardiovascular_Variable = beautify_variable_name(cv_var),
          Dietary_Variable = beautify_group_name(diet_var),
          Sample_Size = ancova_result$n,
          Group1_Name = group_names[1],
          Group1_Mean = round(group_means[1], 3),
          Group1_SE = round(group_se[1], 3),
          Group1_Mean_SE = group1_mean_se,
          Group2_Name = group_names[2],
          Group2_Mean = round(group_means[2], 3),
          Group2_SE = round(group_se[2], 3),
          Group2_Mean_SE = group2_mean_se,
          P_Value = round(ancova_result$overall_p, 4),
          Effect_Size = round(ancova_result$effect_size, 4),
          Significance = significance,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# 保存详细汇总结果（包含均值±SE）
detailed_summary_file <- file.path(output_dir, "detailed_summary_with_means_SE.csv")
write.csv(summary_results, detailed_summary_file, row.names = FALSE)
created_files <- c(created_files, "detailed_summary_with_means_SE.csv")

# 创建简化版汇总表（只保留主要信息）
simple_summary <- summary_results %>%
  select(Cardiovascular_Variable, Dietary_Variable, Sample_Size, 
         Group1_Name, Group1_Mean_SE, Group2_Name, Group2_Mean_SE, 
         P_Value, Significance) %>%
  rename(
    Variable = Cardiovascular_Variable,
    Grouping = Dietary_Variable,
    N = Sample_Size,
    Group_1 = Group1_Name,
    "Group_1_Mean±SE" = Group1_Mean_SE,
    Group_2 = Group2_Name,
    "Group_2_Mean±SE" = Group2_Mean_SE,
    "P_value" = P_Value,
    Sig = Significance
  )

# 保存简化汇总结果
simple_summary_file <- file.path(output_dir, "simple_summary_means_SE.csv")
write.csv(simple_summary, simple_summary_file, row.names = FALSE)
created_files <- c(created_files, "simple_summary_means_SE.csv")

# 按饮食变量分组创建分别的汇总表
for(diet_var in dietary_vars) {
  
  # 筛选该饮食变量的结果
  diet_specific_results <- summary_results %>%
    filter(Dietary_Variable == beautify_group_name(diet_var)) %>%
    select(Cardiovascular_Variable, Sample_Size, 
           Group1_Name, Group1_Mean_SE, Group2_Name, Group2_Mean_SE, 
           P_Value, Effect_Size, Significance) %>%
    arrange(P_Value)
  
  # 保存该饮食变量特定的结果
  diet_file <- file.path(output_dir, paste0("summary_", diet_var, "_means_SE.csv"))
  write.csv(diet_specific_results, diet_file, row.names = FALSE)
  created_files <- c(created_files, basename(diet_file))
  
  cat("✓ 创建", beautify_group_name(diet_var), "专门汇总表:", basename(diet_file), "\n")
}

# 创建显著结果的专门表格
significant_results <- summary_results %>%
  filter(P_Value < 0.05) %>%
  arrange(P_Value) %>%
  select(Cardiovascular_Variable, Dietary_Variable, Sample_Size,
         Group1_Name, Group1_Mean_SE, Group2_Name, Group2_Mean_SE,
         P_Value, Effect_Size, Significance)

if(nrow(significant_results) > 0) {
  sig_file <- file.path(output_dir, "significant_results_means_SE.csv")
  write.csv(significant_results, sig_file, row.names = FALSE)
  created_files <- c(created_files, "significant_results_means_SE.csv")
  cat("✓ 创建显著结果汇总表:", basename(sig_file), "\n")
}

# 输出最终汇总信息
cat("\n=== 分析完成汇总 ===\n")
cat("✓ 总共分析了", nrow(summary_results), "个变量组合\n")
cat("✓ 显著结果 (p<0.05):", sum(summary_results$P_Value < 0.05, na.rm = TRUE), "个\n")
cat("✓ 创建的文件:", length(created_files), "个\n")
cat("✓ 输出目录:", output_dir, "\n")

# 显示最显著的结果（包含均值±SE）
if(nrow(summary_results) > 0) {
  cat("\n前5个最小p值的结果 (包含均值±SE):\n")
  top_results <- summary_results %>%
    arrange(P_Value) %>%
    head(5) %>%
    select(Cardiovascular_Variable, Dietary_Variable, 
           Group1_Name, Group1_Mean_SE, Group2_Name, Group2_Mean_SE, 
           P_Value, Significance)
  
  print(top_results, row.names = FALSE)
}

# 如果有显著结果，单独显示
if(nrow(significant_results) > 0) {
  cat("\n所有显著关联 (p<0.05) 的均值±SE:\n")
  print(significant_results %>%
          select(Cardiovascular_Variable, Dietary_Variable,
                 Group1_Name, Group1_Mean_SE, Group2_Name, Group2_Mean_SE,
                 P_Value, Significance), 
        row.names = FALSE)
}

cat("\n=== 创建的CSV文件说明 ===\n")
cat("1. detailed_summary_with_means_SE.csv - 完整详细结果（包含所有统计信息）\n")
cat("2. simple_summary_means_SE.csv - 简化结果（主要信息）\n")
cat("3. summary_eattime01_means_SE.csv - 进食时长专门结果\n")
cat("4. summary_eatt301_means_SE.csv - 晚餐时间专门结果\n")
if(nrow(significant_results) > 0) {
  cat("5. significant_results_means_SE.csv - 仅显著结果\n")
}

cat("\n=== 分析完成 ===\n")