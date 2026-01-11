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

# 加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)
library(tibble)
library(SingleCellExperiment)
library(Matrix)
library(openxlsx)
library(ggpubr)
library(ggrepel)      # 用于优化标签显示，防止重叠
# devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis) # 用于绘制抖动火山图（jjVolcano）

# 设置工作目录
setwd("文件目录")

#### 不同组别中同一细胞类型的差异表达 (Disease vs Control) ####
#### 数据读取与初始化 ####
# 读取已经经过细胞注释的Seurat对象
sc_data <- readRDS("E:\\scRNA-seq Analysis\\data\\Seurat_celltype.rds")

# 1. 设置细胞类型信息
Idents(sc_data) <- "celltype.3rd"            # 输入细胞注释列的列名，注意修改

# 【选项 A】分析所有亚型 (默认)
type <- unique(sc_data$celltype.3rd)         # 获取所有细胞类型名称列表，用于后续循环

# 【选项 B】只分析特定的几个亚型 (若只想对于某一/多种特定细胞进行差异基因分析)
type <- c("Excitatory", "Inhibitory") 

# 2. 设置分组信息（模糊匹配）
# 定义疾病组的样本ID
disease_pattern <- "GSM9101266|GSM9101267|GSM9101268"

# 新增 group 列，并判断组别
sc_data$group <- ifelse(grepl(disease_pattern, sc_data$orig.ident), "Disease", "Control")

# 设置为因子(Factor)类型数据
sc_data$group <- factor(sc_data$group, levels = c("Control", "Disease"))

# 查看分组结果
print(table(sc_data$orig.ident, sc_data$group))
print(table(sc_data$group))

# 初始化用于存储所有差异分析结果的空数据框
r.deg <- data.frame()

#### 循环计算差异基因 ####
group_case <- "Disease"  # 实验组，注意修改
group_ctrl <- "Control"  # 对照组，注意修改

# 创建结果保存的子文件夹
result_dir <- "DEG_by_Group_results" 
if (!dir.exists(result_dir)) dir.create(result_dir)

# 循环分析
# 注意，此方法需要一个细胞类型中细胞数大于3，否则不会进行差异基因的计算
# 不过如果细胞过少也没必要进行组间对比了，毕竟已经形成新的类型/亚型了
for (i in 1:length(type)) {
  cell_name <- type[i]
  message(paste0("正在分析: ", cell_name, " (", i, "/", length(type), ")"))
  
  try({
    # 1. 计算差异基因
    deg <- FindMarkers(sc_data, 
                       only.pos = FALSE,             # 保留上调和下调
                       min.pct = 0.25,               # 最小表达比例阈值
                       ident.1 = group_case,         # 实验组
                       ident.2 = group_ctrl,         # 对照组
                       group.by = "group",           # 分组列名
                       subset.ident = cell_name,     # 当前循环的细胞类型
                       logfc.threshold = 0.25,       # LogFC筛选阈值
                       test.use = 'wilcox')          # 检验方法
    
    # 2. 在循环内提取基因名
    deg$gene <- rownames(deg)
    
    # 3. 添加注释列
    deg$cluster <- cell_name  # 明确标记细胞类型，供jjVolcano使用
    deg$unm <- i - 1          # 索引标记
    
    # 4. 保存单个细胞类型的中间结果 (文件名去除特殊字符防止报错)
    safe_name <- gsub("/", "_", cell_name)
    file_name <- paste0("0", i, ".", safe_name, '_deg.csv')
    full_path <- file.path(result_dir, file_name)
    write.csv(deg, file = full_path, row.names = FALSE)
    
    # 5. 合并结果
    r.deg <- rbind(r.deg, deg)
  })
}

# 检查合并后的数据结构
head(r.deg)

#### 数据筛选与Top基因提取 ####
# 1. 设定阈值，筛选显著差异基因
s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
# (P.adj < 0.05 且 |Log2FC| > 0.5)

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

#### 绘制抖动火山图 ####
# 计算需要的颜色数量
n_types <- length(unique(r.deg$cluster))

# 预设的基础色板
base_colors <- c("#ec8324","#057c4f","#b7183a","#03957b","#16adc2",
                 "#dd4a3d","#375289","#b5a578","#544f8a",
                 "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")

# 使用 colorRampPalette 进行动态生成
my_tile_col <- colorRampPalette(base_colors)(n_types)

# 绘图
p <- jjVolcano(diffData = r.deg, 
               tile.col = my_tile_col,     # 设置每个cluster的背景色
               pSize = 0.4,                # 基因点的大小
               legend.position = c(0.1, 0.9), 
               celltypeSize = 5,           # X轴细胞类型标签大小
               topGeneN = 0) +             # 设为0，屏蔽默认标签
  geom_text_repel(data = top_label, 
                  aes(x = cluster, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), # 配合抖动散点的位置
                  show.legend = FALSE, 
                  size = 3.5,               # 标签字号
                  box.padding = unit(0.3, "lines"),
                  min.segment.length = 0,   # 总是显示连线
                  max.overlaps = 30)        # 允许更多标签重叠

# 输出与保存
print(p)
ggsave(file.path(result_dir, "Multigroup_Volcano_Plot.pdf"), plot = p, width = 12, height = 8)


