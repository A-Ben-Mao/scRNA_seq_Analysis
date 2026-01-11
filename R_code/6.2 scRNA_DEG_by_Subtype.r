# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面
# 再完成了前面一些列初始的操作后，并完成了细胞注释
# 接下来就可以开始下游的其他分析

# 差异基因分析常用于以下两种情况
# 1.不同组别中同一细胞类型的差异表达
#   如：药物干预 VS 对照、疾病 VS 正常，etc.
# 2.同一细胞类型中不同亚型的差异表达
#   如：新发现的一种亚型有什么特殊的表达谱，etc.

# 1. 加载必要的R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)
library(tibble)
library(ggrepel)
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

# 2. 设置工作目录
setwd("文件目录")

#### 1. 一种Subtype vs. 其他所有Subtypes (One vs Rest) ####
# 3. 数据读取与初始化
# 读取 Seurat 对象
sc_data <- readRDS("E:\\scRNA-seq Analysis\\data\\Seurat_celltype.rds")

# 设置细胞身份 (亚型所在的列名)
Idents(sc_data) <- "celltype.3rd"  # <--- 请修改为你的亚型列名

# 4. 定义要分析的亚型
# 获取数据中所有存在的亚型名称
all_subtypes <- unique(sc_data$celltype.3rd)
all_subtypes

# 【选项 A】分析所有亚型 (默认)
target_types <- all_subtypes

# 【选项 B】只分析特定的几个亚型
target_types <- c("Inhibitory", "Excitatory")

message(paste0("即将分析的亚型数量: ", length(target_types)))

# 创建结果保存目录
result_dir <- "Subtype_OneVsRest_Results"
if (!dir.exists(result_dir)) dir.create(result_dir)

# 初始化用于存储所有结果的大表（用于最后画图）
r.deg <- data.frame()

# 5. 循环计算差异基因 (One vs Rest)
for (i in 1:length(target_types)) {
  cell_name <- target_types[i]
  message(paste0("正在分析 [", i, "/", length(target_types), "]: ", cell_name, " vs Others"))
  
  try({
    # 计算差异基因
    # ident.1 = 当前亚型
    # ident.2 = NULL (默认不写就是 NULL，代表对比剩余所有细胞)
    deg <- FindMarkers(sc_data, 
                       ident.1 = cell_name,
                       ident.2 = NULL,           # 关键：NULL 代表 vs Rest
                       min.pct = 0.25,           # 最小表达比例
                       logfc.threshold = 0.25,   # logFC 阈值
                       test.use = "wilcox",      # 检验方法
                       only.pos = FALSE)         # FALSE 保留负值以便画火山图
    
    # --- 数据整理 ---
    # 1. 提取基因名 (FindMarkers 结果基因名在 rownames)
    deg$gene <- rownames(deg)
    
    # 2. 添加亚型标签 (供绘图分组使用)
    deg$cluster <- cell_name 
    
    # 3. 添加索引 (可选)
    deg$unm <- i - 1
    
    # --- 保存单个结果文件 ---
    # 文件名处理：替换特殊字符防止文件名非法
    safe_name <- gsub("[/|:?]", "_", cell_name)
    file_name <- paste0("OneVsRest_", safe_name, ".csv")
    write.csv(deg, file = file.path(result_dir, file_name), row.names = FALSE)
    
    # --- 合并到总表 ---
    r.deg <- rbind(r.deg, deg)
  })
}

# 检查结果
head(r.deg)

# 数据筛选与标签提取 (为绘图做准备)
# 1. 筛选显著差异基因
# 这里定义显著为: FDR < 0.05 且 |Log2FC| > 0.5
s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

# 2. 添加分类标签
# 标记上下调
s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 0.5 , 'Up', 'Down'))
# 标记显著性等级
s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.05 , 'Highly', 'Lowly'))
# 组合标签
s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)

# 3. 提取绘图用的Top标签 (每个细胞类型取 上调前5 和 下调前5)
top_up_label <- s.deg %>% 
  filter(threshold == "Up") %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()

top_down_label <- s.deg %>% 
  filter(threshold == "Down") %>% 
  group_by(cluster) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% # 负值取最小，即下调最显著
  as.data.frame()

top_label <- rbind(top_up_label, top_down_label)

# 7. 绘制抖动火山图 (jjVolcano)
# 动态生成颜色板
n_types <- length(unique(s.deg$cluster))

# 使用 colorRampPalette 扩展颜色，防止亚型太多颜色不够
my_tile_col <- colorRampPalette(c("#E64B35FF", "#4DBBD5FF", "#00A087FF", 
                                  "#3C5488FF", "#F39B7FFF", "#8491B4FF", 
                                  "#91D1C2FF", "#DC0000FF", "#7E6148FF"))(n_types)

# 绘图
p <- jjVolcano(diffData = r.deg, 
               tile.col = my_tile_col,     # 背景色
               pSize = 0.4,                # 散点大小
               legend.position = c(0.1, 0.9), 
               celltypeSize = 5,           # X轴标签字体大小
               topGeneN = 0) +             # 设置为0，屏蔽默认标签，使用自定义repel
  geom_text_repel(data = top_label, 
                  aes(x = cluster, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), # 配合抖动
                  show.legend = FALSE, 
                  size = 3,                 # 标签字体大小
                  box.padding = unit(0.3, "lines"),
                  min.segment.length = 0,   # 总是显示连线
                  max.overlaps = 50,        # 允许重叠数量
                  fontface = "italic") +    # 基因名斜体
  labs(title = "Subtype Specific Markers (One vs Rest)",
       subtitle = paste0("Comparison: Each Subtype vs All Other Cells"))

# 输出在屏幕
print(p)

# 保存图片
ggsave(file.path(result_dir, "All_Subtypes_OneVsRest_Volcano.pdf"), 
       plot = p, width = 14, height = 8)

#### 2. 指定两个特定亚型直接对比 (One vs One) ####
# 可能用途示例：
# 比如，"Excitatory" 和 "Inhibitory" 的直接差异
# 比如，两种特定亚型细胞的表达差异比较
# 这种对比通常用于画普通火山图，而不是抖动火山图

# 读取Seurat对象
sc_data <- readRDS("E:\\scRNA-seq Analysis\\data\\Seurat_celltype.rds")

# 1. 设置细胞类型信息
Idents(sc_data) <- "celltype.3rd" # 输入细胞注释列的列名，注意修改

# 创建专门保存 One vs One 结果的文件夹
vs_dir <- "Subtype_OneVsOne_Results"
if (!dir.exists(vs_dir)) dir.create(vs_dir)

# 1. 指定对比
ident_1 <- "Excitatory"
ident_2 <- "Inhibitory"

# 差异基因分析
deg_pair <- FindMarkers(sc_data, 
                        ident.1 = ident_1, 
                        ident.2 = ident_2,
                        min.pct = 0.25,
                        logfc.threshold = 0.25)
deg_pair$gene <- rownames(deg_pair)

# 保存 CSV 表格
csv_name <- paste0(ident_1, "_vs_", ident_2, "_diff_genes.csv")
write.csv(deg_pair, file = file.path(vs_dir, csv_name))

# 结果可视化
# 定义显著性阈值
logfc_cutoff <- 0.5
pval_cutoff <- 0.05

# 添加分组标签 (Up, Down, Stable)
deg_pair$direction <- "Stable"
deg_pair$direction[deg_pair$p_val_adj < pval_cutoff & deg_pair$avg_log2FC > logfc_cutoff] <- "Up"
deg_pair$direction[deg_pair$p_val_adj < pval_cutoff & deg_pair$avg_log2FC < -logfc_cutoff] <- "Down"

# 转化为因子，控制绘图层级
deg_pair$direction <- factor(deg_pair$direction, levels = c("Up", "Down", "Stable"))
deg_pair <- deg_pair %>% arrange(desc(direction)) # 排序，让有颜色的点最后画

# 提取 Top 10 标签
# 提取上调前10 和 下调前10 (基于 avg_log2FC)
top_genes <- deg_pair %>%
  filter(direction != "Stable") %>%  # 只在显著基因里找
  group_by(direction) %>%
  slice_max(order_by = abs(avg_log2FC), n = 10) # 上下各取绝对值最大的10个

# 绘制普通火山图
p_volcano <- ggplot(deg_pair, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  # 绘制散点
  geom_point(aes(color = direction), alpha = 0.5, size = 1.5) +
  
  # 自定义颜色 (对应 Up, Down, Stable)
  scale_color_manual(values = c("Up" = "#E64B35FF",      # 红色
                                "Down" = "#00A087FF",    # 绿色/蓝色
                                "Stable" = "grey80")) +  # 灰色
  
  # 添加阈值辅助线
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", size = 0.3) +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", size = 0.3) +
  
  # 添加 Top 10 基因标签
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = 50,
                  segment.color = "black",
                  show.legend = FALSE) +
  
  # 主题设置
  theme_bw() +
  labs(title = paste0(ident_1, " vs ", ident_2),
       subtitle = paste0("Up: ", nrow(subset(deg_pair, direction == "Up")), 
                         " | Down: ", nrow(subset(deg_pair, direction == "Down"))),
       x = "log2 Fold Change", 
       y = "-log10(Adjusted P-value)",
       color = "Significance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# 输出到屏幕
print(p_volcano)

# 7. 保存图片
plot_name <- paste0("Volcano_", ident_1, "_vs_", ident_2, ".pdf")
ggsave(file.path(vs_dir, plot_name), plot = p_volcano, width = 8, height = 7)


