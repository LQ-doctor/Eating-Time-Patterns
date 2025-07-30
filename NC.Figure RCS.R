# 载入必要的包
library(rms)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(eoffice) 
# 确保加载所需的包
library(showtext)
library(grid)

# 根据 Early vs Conv 的比较结果选择差异脂质
rm(list=ls())

# 1. 读取和清理数据
data <- read_excel("TREData.xlsx")
data_clean <- data %>%
  drop_na(eattime, HBP) %>%
  mutate(eattime = case_when(
    eattime < 9 ~ 9,
    eattime > 16 ~ 16,
    TRUE ~ eattime
  ))

# 2. 计算关键分位数
p10 <- round(quantile(data_clean$eattime, 0.10), 1)
p50 <- round(quantile(data_clean$eattime, 0.50), 1)
p90 <- round(quantile(data_clean$eattime, 0.90), 1)
print(paste("P10:", p10, "P50:", p50, "P90:", p90))

# 3. 变量因子化
data_clean$HBP <- factor(data_clean$HBP, levels = c("0", "1"))
data_clean$urban <- as.factor(data_clean$urban)
data_clean$sex <- as.factor(data_clean$sex)
data_clean$MARRIAGE <- as.factor(data_clean$MARRIAGE)
data_clean$FEDU <- as.factor(data_clean$FEDU)
data_clean$MEDU <- as.factor(data_clean$MEDU)
data_clean$FOCUP <- as.factor(data_clean$FOCUP)
data_clean$MOCUP <- as.factor(data_clean$MOCUP)
data_clean$INCOME <- as.factor(data_clean$INCOME)
data_clean$fhistory <- as.factor(data_clean$fhistory)
data_clean$age12 <- as.factor(data_clean$age12)
# 4. 设定数据环境，使用P50作为参考值
dd <- datadist(data_clean)
dd$limits$eattime[2] <- p50  # 设置P50为参考值
options(datadist='dd')

# 5. 拟合RCS模型
fit <- lrm(HBP ~ rcs(eattime, 4) + 
             urban + sex + height_zscore + WHtR + fhistory + 
             MARRIAGE + FEDU + MEDU + FOCUP + MOCUP + age12 +
             INCOME + Energy + total_met,
           data=data_clean)
pred_eattime<- seq(min(data_clean$eattime), max(data_clean$eattime), length.out=500)
# 6. 获取模型结果
an <- anova(fit)
p_values <- round(an[,"P"], 3)
RCS <- Predict(fit, eattime=pred_eattime, fun=exp, ref.zero=TRUE)

rcs_df <- as.data.frame(RCS)

# 7. 创建直方图数据
hist_data <- hist(data_clean$eattime, 
                  breaks=seq(9, 16, length.out=30), 
                  plot=FALSE)
max_count <- max(hist_data$counts)

# 8. 创建ggplot图形
p1_TRE <- ggplot() +
  # 添加直方图
  geom_histogram(data=data_clean, 
                 aes(x=eattime, y=..count../max_count*2),
                 breaks=seq(9, 16, length.out=30),
                 fill='#00AFBB',
                 color='black',
                 size=0.5,
                 alpha=0.7) +
  
  # 添加RCS曲线
  geom_line(data=subset(rcs_df, eattime >= 9 & eattime <= 16), 
            aes(x=eattime, y=yhat),
            color='#E24E67',
            size=1,
            linetype='solid') +
  
  # 添加置信区间
  geom_line(data=subset(rcs_df, eattime >= 9 & eattime <= 16),
            aes(x=eattime, y=lower),
            color='#EA6E80',
            size=1,
            linetype='dashed') +
  geom_line(data=subset(rcs_df, eattime >= 9 & eattime <= 16),
            aes(x=eattime, y=upper),
            color='#EA6E80',
            size=1,
            linetype='dashed') +
  
  # 添加参考线
  geom_hline(yintercept=seq(0, 2, 0.5),
             color="gray",
             linetype='dashed') +
  geom_vline(xintercept=p50,
             color="gray",
             linetype='dashed') +
  
  # 设置坐标轴
  scale_x_continuous(
    name="Eating Duration (hours)",
    limits=c(9, 16),
    breaks=seq(9, 16, 1),
    expand=c(0,0),
    sec.axis=sec_axis(
      ~.,
      name=NULL,
      breaks=c(p10, p50, p90),
      labels=c(paste0("P10"), 
               paste0("P50"),
               paste0("P90"))
    )
  ) +
  scale_y_continuous(
    name="OR (95% CI)",
    breaks=seq(0, 2, 0.5),
    limits=c(0, 2),
    sec.axis=sec_axis(~.*max_count/2,
                      name="Frequency",
                      breaks=c(0, 300, 600, 900, 1200),
                      labels=c("0", "300", "600", "900", "1200"))
  ) +
  
  # 添加p值注释
  annotate("text", x=10, 
           y=1.8,
           label=paste0("P-overall ", 
                        ifelse(p_values[1] < 0.001, "< 0.001", p_values[1]),
                        "\nP-non-linear = ",
                        ifelse(p_values[2] < 0.001, "< 0.001", p_values[2])),
           hjust=0) +
  
  # 设置主题
  theme_classic() +
  theme(
    plot.margin = unit(c(5, 5, 1, 5), "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.line = element_line(color="black"),
    axis.text = element_text(color="black"),
    axis.title = element_text(size=11),
    axis.line.x.top = element_line(color="black"),
    axis.text.x.top = element_text(color="black"),
    axis.ticks.x.top = element_line(color="black"),
    legend.position = "none"
  )

# 显示图形
print(p1_TRE)

topptx(filename = "RCSplot3.pptx",
       width = 6, # 输出图形的宽度和高度
       height = 4.5,append = FALSE,devsize = FALSE,
       units = "in" )# 图形宽度和高度的单位)










# 数据预处理
data <- read_excel("D:/TRE paper/TREData.xlsx")
data_clean <- data %>%
  drop_na(eatt3, HBP) %>%
  mutate(eatt3 = case_when(
    eatt3 < 17 ~ 17,
    eatt3 > 23 ~ 23,
    TRUE ~ eatt3
  ))


# 计算关键分位数
p10 <- round(quantile(data_clean$eatt3, 0.10), 1)
p50 <- round(quantile(data_clean$eatt3, 0.50), 1)
p90 <- round(quantile(data_clean$eatt3, 0.90), 1)
print(paste("P10:", p10, "P50:", p50, "P90:", p90))

# 变量因子化
data_clean$HBP <- factor(data_clean$HBP, levels = c("0", "1"))
data_clean$urban <- as.factor(data_clean$urban)
data_clean$sex <- as.factor(data_clean$sex)
data_clean$MARRIAGE <- as.factor(data_clean$MARRIAGE)
data_clean$FEDU <- as.factor(data_clean$FEDU)
data_clean$MEDU <- as.factor(data_clean$MEDU)
data_clean$FOCUP <- as.factor(data_clean$FOCUP)
data_clean$MOCUP <- as.factor(data_clean$MOCUP)
data_clean$INCOME <- as.factor(data_clean$INCOME)
data_clean$fhistory <- as.factor(data_clean$fhistory)
data_clean$age12 <- as.factor(data_clean$age12)
# 设定数据环境，使用P50作为参考值
dd <- datadist(data_clean)
dd$limits$eatt3[2] <- p50  # 设置P50为参考值
options(datadist='dd')

# 拟合RCS模型
fit <- lrm(HBP ~ rcs(eatt3, 4) + 
             urban + sex + height_zscore + WHtR + fhistory + age12 +
             MARRIAGE + FEDU + MEDU + FOCUP + MOCUP + 
             INCOME + Energy + total_met,
           data=data_clean)

# 获取模型结果
an <- anova(fit)
p_values <- round(an[,"P"], 3)
# 创建更密集的预测点序列
pred_eatt3 <- seq(min(data_clean$eatt3), max(data_clean$eatt3), length.out=1000)
RCS <- Predict(fit, eatt3=pred_eatt3, fun=exp, ref.zero=TRUE)
rcs_df <- as.data.frame(RCS)

# 创建直方图数据
hist_data <- hist(data_clean$eatt3, 
                  breaks=seq(17, 23, length.out=30), 
                  plot=FALSE)
max_count <- max(hist_data$counts)

# 创建ggplot图形
p2_TRE <- ggplot() +
  # 添加直方图
  geom_histogram(data=data_clean, 
                 aes(x=eatt3, y=..count../max_count*2),
                 breaks=seq(17, 23, length.out=30),
                 fill='#00AFBB',
                 color='black',
                 alpha=0.7) +
  
  # 添加RCS曲线
  geom_line(data=subset(rcs_df, eatt3 >= 17 & eatt3 <= 23), 
            aes(x=eatt3, y=yhat),
            color='#E24E67',
            size=1,
            linetype='solid') +
  
  # 添加置信区间
  geom_line(data=subset(rcs_df, eatt3 >= 17 & eatt3 <= 23),
            aes(x=eatt3, y=lower),
            color='#EA6E80',
            size=1,
            linetype='dashed') +
  geom_line(data=subset(rcs_df, eatt3 >= 17 & eatt3 <= 23),
            aes(x=eatt3, y=upper),
            color='#EA6E80',
            size=1,
            linetype='dashed') +
  
  # 添加参考线
  geom_hline(yintercept=seq(0, 2, 0.5),
             color="gray",
             linetype='dashed') +
  geom_vline(xintercept=p50,
             color="gray",
             linetype='dashed') +
  
  # 设置坐标轴
  scale_x_continuous(
    name="Last meal time (hours)",
    limits=c(17, 23),
    breaks=seq(17, 23, 1),
    expand=c(0,0),
    sec.axis=sec_axis(
      ~.,
      name=NULL,
      breaks=c(p10, p50, p90),
      labels=c(paste0("P10"), 
               paste0("P50"),
               paste0("P90"))
    )
  ) +
  scale_y_continuous(
    name="OR (95% CI)",
    breaks=seq(0, 2, 0.5),
    limits=c(0, 2),
    sec.axis=sec_axis(~.*max_count/2,
                      name="Frequency",
                      breaks=c(0, 500, 1000, 1500, 2000),
                      labels=c("0", "500", "1000", "1500", "2000"))
  ) +
  
  # 添加p值注释
  annotate("text", x=18, 
           y=1.8,
           label=paste0("P-overall ", 
                        ifelse(p_values[1] < 0.001, "< 0.001", p_values[1]),
                        "\nP-non-linear = ",
                        ifelse(p_values[2] < 0.001, "< 0.001", p_values[2])),
           hjust=0) +
  
  # 设置主题
  theme_classic() +
  theme(
    plot.margin = unit(c(5, 5, 1, 5), "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.line = element_line(color="black"),
    axis.text = element_text(color="black"),
    axis.title = element_text(size=11),
    axis.line.x.top = element_line(color="black"),
    axis.text.x.top = element_text(color="black"),
    axis.ticks.x.top = element_line(color="black"),
    legend.position = "none"
  )

# 显示图形
print(p2_TRE)
# 载入必要的包
library(patchwork)
library(ggplot2)

# 添加标题标识
p1_TRE <- p1_TRE + 
  ggtitle("A") +
  theme(plot.title = element_text(face = "bold", size = 14))

p2_TRE <- p2_TRE + 
  ggtitle("B") +
  theme(plot.title = element_text(face = "bold", size = 14))

# 组合图形
combined_plot <- p1_TRE + p2_TRE +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(theme = theme(plot.margin = unit(c(1, 1, 1, 1), "cm")))

# 显示组合后的图形
print(combined_plot)
# 保存图形
ggsave("RCS_TRE_eatt3.pdf", p1_TRE, width=10, height=8, dpi=300)
ggsave("RCS_TRE_eatt3.tiff", p1_TRE, width=10, height=8, dpi=300)
