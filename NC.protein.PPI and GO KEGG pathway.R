# 读取并处理数据
ppi_df <- read_delim(
  file = "string_interactions.tsv",  # 更正的文件路径
  delim = "\t",
  col_names = TRUE
)

# 查看数据结构
str(ppi_df)
head(ppi_df)

# 创建边列表
edge_list <- ppi_df %>%
  select(node1, node2, combined_score) %>%
  rename(
    from = node1,
    to = node2,
    value = combined_score
  )

# 数据过滤 - 可选，如果网络太大可以过滤掉低分值的互作用
edge_list_filtered <- edge_list %>%
  filter(value > 700)  # 使用700分作为阈值，可以根据需要调整

# 创建图形并计算中心度
graph_df <- as_tbl_graph(edge_list_filtered) %>%
  tidygraph::mutate(
    Popularity = centrality_degree(mode = 'all'),
    Betweenness = centrality_betweenness(),
    Degree_group = case_when(
      Popularity > 20 ~ "High",
      Popularity > 10 ~ "Medium",
      TRUE ~ "Low"
    )
  )

# 查看节点数量
node_count <- graph_df %>% activate(nodes) %>% as_tibble() %>% nrow()
edge_count <- graph_df %>% activate(edges) %>% as_tibble() %>% nrow()
cat("节点数量:", node_count, "\n")
cat("边数量:", edge_count, "\n")

# 绘制网络图
set.seed(123)  # 设置随机种子以确保结果可复现

# 为中等大小的网络使用Fruchterman-Reingold布局
if(node_count < 200) {
  p1 <- ggraph(graph_df, layout = 'fr') +
    geom_edge_link(aes(color = value, width = value/1000), 
                   alpha = 0.7) +
    geom_node_point(aes(size = Popularity,
                        fill = Degree_group), 
                    color = "#000000",
                    alpha = 1, 
                    shape = 21) +
    geom_node_text(aes(label = name), 
                   repel = TRUE, 
                   size = 3,
                   check_overlap = TRUE) +
    scale_edge_colour_distiller(palette = "RdPu", 
                                direction = 1,
                                name = "Combined score") +
    scale_edge_width(range = c(0.1, 1.5), name = "Edge Weight") +
    scale_size_continuous(
      range = c(2, 10),
      name = "Degree") +
    scale_fill_manual(
      values = c("High" = "#d73027", "Medium" = "#fc8d59", "Low" = "#fee090"),
      name = "Connectivity"
    ) +
    theme_graph() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(title = "Protein-Protein Interaction Network")
} else {
  # 为大型网络使用环形布局
  p1 <- ggraph(graph_df, layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(color = value, width = value/1000), 
                  alpha = 0.5) +
    geom_node_point(aes(size = Popularity,
                        fill = Degree_group), 
                    color = "#000000",
                    alpha = 1, 
                    shape = 21) +
    geom_node_text(
      aes(label = ifelse(Popularity > quantile(Popularity, 0.9), name, "")),
      repel = TRUE,
      size = 3
    ) +
    scale_edge_colour_distiller(palette = "RdPu", 
                                direction = 1,
                                name = "Combined score") +
    scale_edge_width(range = c(0.1, 1), name = "Edge Weight") +
    scale_size_continuous(
      range = c(1, 8),
      name = "Degree") +
    scale_fill_manual(
      values = c("High" = "#d73027", "Medium" = "#fc8d59", "Low" = "#fee090"),
      name = "Connectivity"
    ) +
    theme_graph() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(title = "Protein-Protein Interaction Network") +
    coord_fixed()
}

# 显示网络图
print(p1)

# 保存网络图
ggsave("ppi_network.pdf", p1, width = 12, height = 10)
ggsave("ppi_network.png", p1, width = 12, height = 10, dpi = 300)

# 导出到PPT
topptx(p1,
       filename = "PPI_Network.pptx",
       width = 12, 
       height = 10,
       append = FALSE,
       devsize = FALSE,
       units = "in")

# 找出中心度最高的前20个蛋白质
top_20_proteins <- graph_df %>%
  activate(nodes) %>%
  as_tibble() %>%
  arrange(desc(Popularity)) %>%
  slice_head(n = 20)

# 显示中心度最高的蛋白质
print("Top 20 proteins by connectivity degree:")
print(top_20_proteins %>% 
        select(name, Popularity, Betweenness, Degree_group) %>%
        arrange(desc(Popularity)))

# 保存中心度最高的蛋白质列表
write.csv(top_20_proteins, "top_20_hub_proteins.csv", row.names = FALSE)

# 可选：为高中心度蛋白质创建子网络
if(node_count > 50) {
  # 获取高中心度蛋白质名称
  hub_proteins <- top_20_proteins$name
  
  # 创建子网络
  subgraph <- graph_df %>%
    convert(to_subgraph, name %in% hub_proteins | 
              .N()$name[from] %in% hub_proteins & .N()$name[to] %in% hub_proteins)
  
  # 绘制子网络
  p2 <- ggraph(subgraph, layout = 'fr') +
    geom_edge_link(aes(color = value, width = value/1000), 
                   alpha = 0.7) +
    geom_node_point(aes(size = Popularity,
                        fill = Degree_group), 
                    color = "#000000",
                    alpha = 1, 
                    shape = 21) +
    geom_node_text(aes(label = name), 
                   repel = TRUE, 
                   size = 3.5) +
    scale_edge_colour_distiller(palette = "RdPu", 
                                direction = 1,
                                name = "Combined score") +
    scale_edge_width(range = c(0.2, 2), name = "Edge Weight") +
    scale_size_continuous(
      range = c(3, 12),
      name = "Degree") +
    scale_fill_manual(
      values = c("High" = "#d73027", "Medium" = "#fc8d59", "Low" = "#fee090"),
      name = "Connectivity"
    ) +
    theme_graph() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(title = "Hub Proteins Interaction Network")
  
  # 显示子网络
  print(p2)
  
  # 保存子网络图
  ggsave("hub_ppi_network.pdf", p2, width = 10, height = 8)
  ggsave("hub_ppi_network.png", p2, width = 10, height = 8, dpi = 300)
  
  # 导出到PPT
  topptx(p2,
         filename = "Hub_PPI_Network.pptx",
         width = 10,
         height = 8,
         append = FALSE,
         devsize = FALSE,
         units = "in")
}





# 加载必要的包
library(tidyverse)
library(ggplot2)
library(ggnewscale)  # 用于多重图例
library(openxlsx)    # 用于读取Excel文件

# 读取数据
enrichment_data <- read.xlsx("D:/TREpaper/通路富集.xlsx")

# 处理数据
df <- enrichment_data %>%
  # 提取通路类型和比较组
  mutate(
    ONTOLOGY = case_when(
      grepl("BP", Category) ~ "BP",
      grepl("CC", Category) ~ "CC",
      grepl("MF", Category) ~ "MF",
      TRUE ~ as.character(Category)
    ),
    Group = comparison,  # 保持原始比较组标签
    GeneRatio = Generatio  # 使用Generatio列
  ) %>%
  # 移除KEGG数据
  filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 获取所有比较组
all_groups <- unique(df$Group)

# 计算每个Term在每个比较组中是否存在
term_presence <- df %>%
  select(ONTOLOGY, Term, Group, GeneRatio, PValue) %>%
  distinct() %>%
  group_by(ONTOLOGY, Term) %>%
  summarise(
    groups_present = list(unique(Group)),
    group_count = n_distinct(Group),
    all_groups_present = n_distinct(Group) == length(all_groups),
    mean_generatio = mean(GeneRatio),
    min_pvalue = min(PValue),
    .groups = 'drop'
  ) %>%
  arrange(ONTOLOGY, desc(group_count), desc(mean_generatio), min_pvalue)

# 第一步：选择在所有比较组中都存在的Term
terms_in_all_groups <- term_presence %>%
  filter(all_groups_present == TRUE) %>%
  select(ONTOLOGY, Term)

# 第二步：对于每个ONTOLOGY，如果在所有组中存在的Term少于10个，则补充
terms_to_supplement <- data.frame()

for(o in unique(df$ONTOLOGY)) {
  # 获取该ONTOLOGY下在所有比较组中都存在的Term数量
  all_groups_count <- term_presence %>%
    filter(ONTOLOGY == o, all_groups_present == TRUE) %>%
    nrow()
  
  # 如果少于10个，则需要补充
  if(all_groups_count < 10) {
    # 需要补充的数量
    to_add_count <- 10 - all_groups_count
    
    # 已选的Term
    selected_terms <- term_presence %>%
      filter(ONTOLOGY == o, all_groups_present == TRUE) %>%
      pull(Term)
    
    # 按照组合排序选择额外的Term (先按出现组数、再按GeneRatio、最后按P值)
    additional_terms <- term_presence %>%
      filter(ONTOLOGY == o, !Term %in% selected_terms) %>%
      arrange(desc(group_count), desc(mean_generatio), min_pvalue) %>%
      head(to_add_count)
    
    # 合并到要补充的列表中
    terms_to_supplement <- bind_rows(
      terms_to_supplement,
      additional_terms %>% select(ONTOLOGY, Term)
    )
  }
}

# 合并所有要展示的Term
all_selected_terms <- bind_rows(
  terms_in_all_groups %>% mutate(source = "all_groups"),
  terms_to_supplement %>% mutate(source = "supplemental")
)

# 从原始数据中筛选出这些selected terms的数据
df_filtered <- df %>%
  inner_join(all_selected_terms, by = c("ONTOLOGY", "Term"))

# 处理Term名称
df_filtered <- df_filtered %>%
  mutate(
    Description = case_when(
      grepl("~", Term) ~ str_extract(Term, "(?<=~).*$"),
      grepl("`", Term) ~ str_remove(str_extract(Term, "`.*$"), "`"),
      TRUE ~ Term
    )
  )

# 创建Term排序 (先按source优先展示所有组都有的，然后按GeneRatio和P值)
term_order <- df_filtered %>%
  group_by(ONTOLOGY, Term, source) %>%
  summarise(
    mean_generatio = mean(GeneRatio),
    min_pvalue = min(PValue),
    .groups = "drop"
  ) %>%
  mutate(source_priority = ifelse(source == "all_groups", 1, 0)) %>%
  arrange(ONTOLOGY, desc(source_priority), desc(mean_generatio), min_pvalue) %>%
  left_join(df_filtered %>% select(ONTOLOGY, Term, Description) %>% distinct(), 
            by = c("ONTOLOGY", "Term"))

# 创建最终的绘图数据框
df_plot <- df_filtered %>%
  mutate(
    Group = factor(Group, levels = c("Early_TRE_vs_Extend_EW", "Late_TRE_vs_Extend_EW", "Early_TRE_vs_Late_TRE")),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF")),
    Description = factor(Description, levels = unique(term_order$Description))
  )

# 创建各组数据框
df1 <- df_plot %>% filter(Group == "Early_TRE_vs_Extend_EW")
df2 <- df_plot %>% filter(Group == "Late_TRE_vs_Extend_EW")
df3 <- df_plot %>% filter(Group == "Early_TRE_vs_Late_TRE")

# 绘图
p <- ggplot(data = df_plot, aes(x = Group, y = Description)) +
  # Early_TRE_vs_Extend_EW
  geom_point(data = df1,
             aes(size = GeneRatio, color = -log10(PValue)),
             shape = 16,
             alpha = 0.8) +
  scale_color_gradient(
    low = "#d4b9da", 
    high = "#ce1256",
    name = "-log10(P-value)\nEarly TRE vs EW",
    na.value = NA
  ) +
  
  new_scale_color() +
  # Late_TRE_vs_Extend_EW
  geom_point(data = df2,
             aes(size = GeneRatio, color = -log10(PValue)),
             shape = 18,
             alpha = 0.8) +
  scale_color_gradient(
    low = "#c7e9c0", 
    high = "#006d2c",
    name = "-log10(P-value)\nLate TRE vs EW",
    na.value = NA
  ) +
  
  new_scale_color() +
  # Early_TRE_vs_Late_TRE
  geom_point(data = df3,
             aes(size = GeneRatio, color = -log10(PValue)),
             shape = 15,
             alpha = 0.8) +
  scale_color_gradient(
    low = "#a6cee3", 
    high = "#2b8cbe",
    name = "-log10(P-value)\nEarly TRE vs Late TRE",
    na.value = NA
  ) +
  
  # 设置大小比例尺
  scale_size(name = "Gene Ratio", 
             range = c(2, 6),
             guide = guide_legend(override.aes = list(alpha = 1))) +
  
  # 添加分面
  facet_grid(ONTOLOGY ~ ., 
             scales = "free_y",
             space = "free_y") +
  
  scale_y_discrete(position = "right") +
  
  labs(x = "", y = "",
       title = "GO Term Enrichment Analysis") +
  
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white"),
    
    legend.box = "vertical",
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.8, "lines"),
    
    plot.margin = ggplot2::margin(t = 15, r = 30, b = 15, l = 15, unit = "pt")
  )

# 显示图形
print(p)

# 保存图形到文件
# 使用eoffice保存到PPT
topptx(p, 
       filename = "pro_Results/GO_enrichment_plot.pptx",
       width = 6, 
       height = 8,
       append = FALSE,
       devsize = FALSE,
       units = "in")



# 设置工作目录
setwd("D:/TREpaper")  # 替换为您的工作目录

# 创建结果目录 
dir.create("pro_Results", showWarnings = FALSE)

# 加载必要的包
library(tidyverse)
library(ggplot2)
library(ggnewscale)  # 用于多重图例
library(openxlsx)    # 用于读取Excel文件

# 读取数据
enrichment_data <- read.xlsx("D:/TREpaper/通路富集.xlsx")

# 处理数据
df <- enrichment_data %>%
  # 提取通路类型和比较组
  mutate(
    ONTOLOGY = case_when(
      grepl("BP", Category) ~ "BP",
      grepl("CC", Category) ~ "CC",
      grepl("MF", Category) ~ "MF",
      TRUE ~ as.character(Category)
    ),
    Group = comparison,  # 保持原始比较组标签
    GeneRatio = Generatio  # 使用Generatio列
  ) %>%
  # 移除KEGG数据
  filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 获取所有比较组
all_groups <- unique(df$Group)

# 计算每个Term在每个比较组中是否存在
term_presence <- df %>%
  select(ONTOLOGY, Term, Group, GeneRatio, PValue) %>%
  distinct() %>%
  group_by(ONTOLOGY, Term) %>%
  summarise(
    groups_present = list(unique(Group)),
    group_count = n_distinct(Group),
    all_groups_present = n_distinct(Group) == length(all_groups),
    mean_generatio = mean(GeneRatio),
    min_pvalue = min(PValue),
    .groups = 'drop'
  ) %>%
  arrange(ONTOLOGY, desc(group_count), desc(mean_generatio), min_pvalue)

# 第一步：选择在所有比较组中都存在的Term
terms_in_all_groups <- term_presence %>%
  filter(all_groups_present == TRUE) %>%
  select(ONTOLOGY, Term)

# 第二步：对于每个ONTOLOGY，如果在所有组中存在的Term少于10个，则补充
terms_to_supplement <- data.frame()

for(o in unique(df$ONTOLOGY)) {
  # 获取该ONTOLOGY下在所有比较组中都存在的Term数量
  all_groups_count <- term_presence %>%
    filter(ONTOLOGY == o, all_groups_present == TRUE) %>%
    nrow()
  
  # 如果少于10个，则需要补充
  if(all_groups_count < 10) {
    # 需要补充的数量
    to_add_count <- 10 - all_groups_count
    
    # 已选的Term
    selected_terms <- term_presence %>%
      filter(ONTOLOGY == o, all_groups_present == TRUE) %>%
      pull(Term)
    
    # 按照组合排序选择额外的Term (先按出现组数、再按GeneRatio、最后按P值)
    additional_terms <- term_presence %>%
      filter(ONTOLOGY == o, !Term %in% selected_terms) %>%
      arrange(desc(group_count), desc(mean_generatio), min_pvalue) %>%
      head(to_add_count)
    
    # 合并到要补充的列表中
    terms_to_supplement <- bind_rows(
      terms_to_supplement,
      additional_terms %>% select(ONTOLOGY, Term)
    )
  }
}

# 合并所有要展示的Term
all_selected_terms <- bind_rows(
  terms_in_all_groups %>% mutate(source = "all_groups"),
  terms_to_supplement %>% mutate(source = "supplemental")
)

# 从原始数据中筛选出这些selected terms的数据
df_filtered <- df %>%
  inner_join(all_selected_terms, by = c("ONTOLOGY", "Term"))

# 处理Term名称
df_filtered <- df_filtered %>%
  mutate(
    Description = case_when(
      grepl("~", Term) ~ str_extract(Term, "(?<=~).*$"),
      grepl("`", Term) ~ str_remove(str_extract(Term, "`.*$"), "`"),
      TRUE ~ Term
    )
  )

# 创建Term排序 (先按source优先展示所有组都有的，然后按GeneRatio和P值)
term_order <- df_filtered %>%
  group_by(ONTOLOGY, Term, source) %>%
  summarise(
    mean_generatio = mean(GeneRatio),
    min_pvalue = min(PValue),
    .groups = "drop"
  ) %>%
  mutate(source_priority = ifelse(source == "all_groups", 1, 0)) %>%
  arrange(ONTOLOGY, desc(source_priority), desc(mean_generatio), min_pvalue) %>%
  left_join(df_filtered %>% select(ONTOLOGY, Term, Description) %>% distinct(), 
            by = c("ONTOLOGY", "Term"))

# 创建最终的绘图数据框
df_plot <- df_filtered %>%
  mutate(
    Group = factor(Group, levels = c("Early_TRE_vs_Extend_EW", "Late_TRE_vs_Extend_EW", "Early_TRE_vs_Late_TRE")),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF")),
    Description = factor(Description, levels = unique(term_order$Description))
  )

# 创建各组数据框
df1 <- df_plot %>% filter(Group == "Early_TRE_vs_Extend_EW")
df2 <- df_plot %>% filter(Group == "Late_TRE_vs_Extend_EW")
df3 <- df_plot %>% filter(Group == "Early_TRE_vs_Late_TRE")

# 绘图
p <- ggplot(data = df_plot, aes(x = Group, y = Description)) +
  # Early_TRE_vs_Extend_EW
  geom_point(data = df1,
             aes(size = GeneRatio, color = -log10(PValue)),
             shape = 16,
             alpha = 0.8) +
  scale_color_gradient(
    low = "#d4b9da", 
    high = "#ce1256",
    name = "-log10(P-value)\nEarly TRE vs EW",
    na.value = NA
  ) +
  
  new_scale_color() +
  # Late_TRE_vs_Extend_EW
  geom_point(data = df2,
             aes(size = GeneRatio, color = -log10(PValue)),
             shape = 18,
             alpha = 0.8) +
  scale_color_gradient(
    low = "#c7e9c0", 
    high = "#006d2c",
    name = "-log10(P-value)\nLate TRE vs EW",
    na.value = NA
  ) +
  
  new_scale_color() +
  # Early_TRE_vs_Late_TRE
  geom_point(data = df3,
             aes(size = GeneRatio, color = -log10(PValue)),
             shape = 15,
             alpha = 0.8) +
  scale_color_gradient(
    low = "#a6cee3", 
    high = "#2b8cbe",
    name = "-log10(P-value)\nEarly TRE vs Late TRE",
    na.value = NA
  ) +
  
  # 设置大小比例尺
  scale_size(name = "Gene Ratio", 
             range = c(2, 6),
             guide = guide_legend(override.aes = list(alpha = 1))) +
  
  # 添加分面
  facet_grid(ONTOLOGY ~ ., 
             scales = "free_y",
             space = "free_y") +
  
  scale_y_discrete(position = "right") +
  
  labs(x = "", y = "",
       title = "GO Term Enrichment Analysis") +
  
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white"),
    
    legend.box = "vertical",
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.8, "lines"),
    
    plot.margin = ggplot2::margin(t = 15, r = 30, b = 15, l = 15, unit = "pt")
  )

# 显示图形
print(p)

# 保存图形到文件
ggsave("pro_Results/GO_enrichment_plot.pdf", p, width = 10, height = 10, units = "in")
ggsave("pro_Results/GO_enrichment_plot.png", p, width = 10, height = 10, units = "in", dpi = 300)
