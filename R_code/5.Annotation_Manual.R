# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao
# 在这里抛砖引玉，希望大家能集万家之所长

#### 写在前面 ####
# 参考同领域专业文献进行注释
# 如果你是根据已经发表文献的二次分析
# 那么可以直接参考原始文献对的细胞类型marker进行注释
# 即使整合了多个样本的二次分析也可参考
# 一般在method、supplement table中可见

# 一些前辈们整理的细胞类型marker
# https://www.jianshu.com/p/15dddefc7038
# https://www.jianshu.com/p/3b28efe61a8f

# 需要注意的是一个Gene可能会有多个名字
# 如免疫细胞的marker：CD45，其在矩阵中的名字为PTPRC
# 可根据GeneCards网站查询：https://www.genecards.org/

#### 主要原则：分级注释 ####
# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(dplyr)

# 设置工作目录
setwd("文件目录")

# 读取已经经过降维的Seurat对象
sc_data = readRDS("Seurat_DimensionalityReduction_harmony_40_0.5.rds")

# 查看各个cluster的细胞数
table(sc_data@meta.data$seurat_clusters)

#### 第一级注释：大类注释 ####
# 注释需要结合测序组织，本次示例的数据是脑组织
major_markers <- list(
  # 1. 免疫细胞 (Immune) - 来源于血液或小胶质细胞
  Immune = c("PTPRC"), # CD45, 所有白细胞
  # 2. 神经/胶质总谱系 (Neuro/Glial)
  Neural_Lineage = c("MAP2", "RBFOX3", "GFAP"), # MAP2/NeuN(神经元), GFAP(胶质)
  # 3. 内皮/血管系统 (Endothelial/Vascular)
  Endothelial = c("PECAM1", "CLDN5", "VWF"), # CD31, 血脑屏障紧密连接蛋白
  # 4. 间质/成纤维细胞 (Stromal/Fibroblasts) - 通常来自脑膜或血管周细胞
  Stromal = c("COL1A1", "PDGFRB", "DCN"), 
  # 5. 上皮细胞 (Epithelial) - 脑中少见，可能是脉络丛上皮
  Epithelial = c("EPCAM", "KRT8", "KRT18"),
  # 6. 增殖细胞 (Proliferating) - 可能是神经前体或肿瘤
  Cycle = c("MKI67", "TOP2A")
)

# 展平列表并去重
major_markers_flat <- unique(unlist(major_markers))
# 检查基因是否存在于数据中
valid_major <- major_markers_flat[major_markers_flat %in% rownames(sc_data)]

# 绘制第一步气泡图
p_major <- DotPlot(sc_data, features = valid_major, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "#FFE4B5", "#FF4500", "#8B0000")) + # 白-黄-橙-深红
  labs(title = "Step 1: Broad Lineage Markers", y = "Clusters", x = "Markers")
print(p_major)

# 绘制第一步cluster点图
pp <- FeaturePlot(sc_data, features = valid_major, ncol = 3)
pp

#### 第一次注释(空白) ####
sc_data$celltype.1st <- recode(sc_data@meta.data$seurat_clusters,
                               "0" = "0",
                               "1" = "1",
                               "2" = "2",
                               "3" = "3",
                               "4" = "4",
                               "5" = "5",
                               "6" = "6",
                               "7" = "7",
                               "8" = "8",
                               "9" = "9",
                               "10" = "10",
                               "11" = "11",
                               "12" = "12",
                               "13" = "13",
                               "14" = "14",
                               "15" = "15",
                               "16" = "16",
                               "17" = "17",
                               "18" = "18",
                               "19" = "19",
                               "20" = "20",
                               "21" = "21",
                               "22" = "22",
                               "23" = "23",
                               "24" = "24",
                               "25" = "25",
                               "26" = "26",
                               "27" = "27",
)

#### 第一次注释(已完成) ####
sc_data$celltype.1st <- recode(sc_data@meta.data$seurat_clusters,
                               "0" = "0",
                               "1" = "1",
                               "2" = "2",
                               "3" = "Immune",
                               "4" = "Neural_Lineage",
                               "5" = "Neural_Lineage",
                               "6" = "Neural_Lineage",
                               "7" = "Neural_Lineage",
                               "8" = "Neural_Lineage",
                               "9" = "Neural_Lineage",
                               "10" = "Neural_Lineage",
                               "11" = "Neural_Lineage",
                               "12" = "Neural_Lineage",
                               "13" = "Immune",
                               "14" = "Immune",
                               "15" = "Neural_Lineage",
                               "16" = "Neural_Lineage",
                               "17" = "Neural_Lineage",
                               "18" = "Neural_Lineage",
                               "19" = "Neural_Lineage",
                               "20" = "Stromal",
                               "21" = "Endothelial",
                               "22" = "Neural_Lineage",
                               "23" = "Neural_Lineage",
                               "24" = "Immune/Neural_Lineage",
                               "25" = "Neural_Lineage",
                               "26" = "Immune",
                               "27" = "Immune/Neural_Lineage",
)

# 查看各个细胞类型的细胞数
table(sc_data@meta.data$celltype.1st)

# 可视化第一次注释的结果
p = DimPlot(sc_data, reduction = "umap", label = T,group.by = "seurat_clusters")
p1 = DimPlot(sc_data, reduction = "umap", label = T,group.by = "celltype.1st")
p|p1

#### 第二级注释：亚群注释 ####
sub_markers <- list(
  # --- 胶质细胞 (Glial Cells) ---
  # 星形胶质细胞 (Astrocytes)
  Astrocytes = c("AQP4", "GFAP", "ALDH1L1", "SLC1A2"),
  # 少突胶质细胞 (Oligodendrocytes) - 形成髓鞘
  Oligodendrocytes = c("MOG", "MBP", "PLP1"),
  # 少突胶质前体细胞 (OPCs)
  OPC = c("PDGFRA", "CSPG4", "OLIG1"), # CSPG4 is NG2
  # 小胶质细胞 (Microglia) - 脑内的免疫细胞
  Microglia = c("AIF1", "CX3CR1", "C1QA", "HEXB", "CSF1R"), # AIF1 is Iba1
  # --- 神经元 (Neurons) ---
  # 神经元通用
  Pan_Neuronal = c("RBFOX3", "SYT1", "SNAP25"), # RBFOX3 is NeuN
  # 兴奋性神经元 (Excitatory / Glutamatergic)
  Excitatory = c("SLC17A7", "CAMK2A", "SATB2"), # SLC17A7 is VGLUT1
  # 抑制性神经元 (Inhibitory / GABAergic)
  Inhibitory = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST"), # SLC32A1 is VGAT
  # --- 其他结构细胞 (Others) ---
  # 室管膜细胞 (Ependymal cells) - 位于脑室
  Ependymal = c("FOXJ1", "DNAH5"),
  # 血管周细胞 (Pericytes)
  Pericytes = c("PDGFRB", "ACTA2", "RGS5")
)

# 展平列表并去重
sub_markers_flat <- unique(unlist(sub_markers))
# 检查基因是否存在于数据中
valid_sub <- sub_markers_flat[sub_markers_flat %in% rownames(sc_data)]

# 绘制第二步气泡图
p_sub <- DotPlot(sc_data, features = valid_sub, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "#E6E6FA", "#483D8B", "#191970")) + # 白-浅紫-深蓝
  labs(title = "Step 2: Brain Specific Subtypes", y = "Clusters", x = "Markers") +
  theme(axis.text.x = element_text(size = 9)) # 字体稍微调小一点以免重叠
print(p_sub)

#### 第二次注释(已完成) ####
sc_data$celltype.2nd <- recode(sc_data@meta.data$seurat_clusters,
                               "0" = "Oligodendrocytes",
                               "1" = "Oligodendrocytes",
                               "2" = "Astrocytes",
                               "3" = "Immune/Microglia",
                               "4" = "Neural_Lineage",
                               "5" = "Neural_Lineage/Excitatory",
                               "6" = "Neural_Lineage/Inhibitory",
                               "7" = "Neural_Lineage/Inhibitory",
                               "8" = "Neural_Lineage/Excitatory",
                               "9" = "Neural_Lineage/Excitatory",
                               "10" = "Neural_Lineage/Inhibitory",
                               "11" = "Neural_Lineage/Excitatory",
                               "12" = "Neural_Lineage/Astrocytes",
                               "13" = "Immune/Microglia",
                               "14" = "Immune",
                               "15" = "Neural_Lineage/Excitatory",
                               "16" = "Neural_Lineage/Inhibitory",
                               "17" = "Neural_Lineage/Oligodendrocytes",
                               "18" = "Neural_Lineage/Excitatory",
                               "19" = "Neural_Lineage/Excitatory",
                               "20" = "Stromal/Pericytes",
                               "21" = "Endothelial",
                               "22" = "Neural_Lineage/Excitatory",
                               "23" = "Neural_Lineage/Inhibitory",
                               "24" = "Immune/Neural_Lineage/Astrocytes",
                               "25" = "Neural_Lineage/Excitatory",
                               "26" = "Immune/Oligodendrocytes/Microglia",
                               "27" = "Immune/Neural_Lineage/Microglia",
)

#### 第二次注释(根据生物学知识简化) ####
sc_data$celltype.2nd <- recode(sc_data@meta.data$seurat_clusters,
                               "0" = "Oligodendrocytes",
                               "1" = "Oligodendrocytes",
                               "2" = "Astrocytes",
                               "3" = "Microglia",
                               "4" = "Neural_Lineage",
                               "5" = "Excitatory",
                               "6" = "Inhibitory",
                               "7" = "Inhibitory",
                               "8" = "Excitatory",
                               "9" = "Excitatory",
                               "10" = "Inhibitory",
                               "11" = "Excitatory",
                               "12" = "Astrocytes",
                               "13" = "Microglia",
                               "14" = "Immune",
                               "15" = "Excitatory",
                               "16" = "Inhibitory",
                               "17" = "Oligodendrocytes",
                               "18" = "Excitatory",
                               "19" = "Excitatory",
                               "20" = "Pericytes",
                               "21" = "Endothelial",
                               "22" = "Excitatory",
                               "23" = "Inhibitory",
                               "24" = "Immune/Neural_Lineage/Astrocytes",
                               "25" = "Excitatory",
                               "26" = "Immune/Oligodendrocytes/Microglia",
                               "27" = "Immune/Neural_Lineage/Microglia",
)

# 查看各个细胞类型的细胞数
table(sc_data@meta.data$celltype.2nd)

# 可视化第一次注释的结果
p2 = DimPlot(sc_data, reduction = "umap", label = T,group.by = "celltype.2nd")
p|p1|p2
# p2

#### FindAllMarkers函数 ####
# Seurat包提供了两种主要的差异表达分析函数：
# FindMarkers()：用于两组细胞间的定向比较
# FindAllMarkers()：自动对所有聚类分组进行差异分析
# 运行时间和细胞数、cluster数呈正相关，整体而言时间较长

# 该的R包可加速wilcox方法，可选择下载
# devtools::install_github('immunogenomics/presto')

# 作用：寻找当前assay下各个identity的marker基因
DefaultAssay(sc_data)
# DefaultAssay(sc_data) = "SCT"
table(Idents(sc_data))
# Idents(sc_data) <- "seurat_clusters" 

?FindAllMarkers
sc_data.markers <- FindAllMarkers(sc_data,
                                  only.pos = TRUE,        # 只保留上调基因
                                  min.pct = 0.25,         # 至少在25%细胞表达
                                  logfc.threshold = 0.25, # 设定logfc阈值
                                  test.use = "wilcox"     # 默认检验方法
)

write.csv(sc_data.markers,file="Markers_logfc0.25.csv")

# 读取已经处理好的Marker文件
sc_data.markers <- read.csv("Markers_logfc0.25.csv",header = TRUE)
rownames(sc_data.markers) <- sc_data.markers[,1]
sc_data.markers = sc_data.markers[,-1]

# 按照log2FC筛选每个cluster的top基因
top10 <- sc_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# 展示top基因的热图：
DoHeatmap(sc_data, features = top10$gene) + NoLegend()

# 转换top基因的表达形式
top10_formatted <- top10 %>%
  group_by(cluster) %>%
  summarise(gene_list = paste(gene, collapse = ",")) %>%
  ungroup()
write.csv(top10_formatted,file="Markers_top.csv")

#### “邪修”技巧 ####
# 技巧1.使用AI辅助注释，提示词参考：
# 我现在正在进行XX组织的单细胞分析，经过降维聚类共有XX簇，
# 以上每一行都是一个细胞簇的top基因，请你辅助我进行单细胞的细胞类型注释。
# 结果以表格的形式呈现，需要包含细胞簇序号、核心Marker基因、推测细胞类型、依据说明、置信度

# 技巧2.使用CellMarker网站辅助注释
# http://www.bio-bigdata.center/

#### 第三次注释(结合AI进行修改) ####
sc_data$celltype.3rd <- recode(sc_data@meta.data$seurat_clusters,
                               "0" = "Oligodendrocytes",
                               "1" = "Oligodendrocytes",
                               "2" = "Astrocytes",
                               "3" = "Microglia",
                               "4" = "Neural_Lineage",
                               "5" = "Excitatory",
                               "6" = "Inhibitory",
                               "7" = "Inhibitory",
                               "8" = "Excitatory",
                               "9" = "Excitatory",
                               "10" = "Inhibitory",
                               "11" = "Excitatory",
                               "12" = "Astrocytes",
                               "13" = "Microglia",
                               "14" = "Immune",
                               "15" = "Excitatory",
                               "16" = "Inhibitory",
                               "17" = "Oligodendrocytes",
                               "18" = "Excitatory",
                               "19" = "Excitatory",
                               "20" = "Pericytes",
                               "21" = "Endothelial",
                               "22" = "Excitatory",
                               "23" = "Inhibitory",
                               "24" = "Immune/Neural_Lineage/Astrocytes",
                               "25" = "Excitatory",
                               "26" = "Immune/Oligodendrocytes/Microglia",
                               "27" = "Immune/Neural_Lineage/Microglia",
)

# 查看各个细胞类型的细胞数
table(sc_data@meta.data$celltype.3rd)

# 可视化第一次注释的结果
p3 = DimPlot(sc_data, reduction = "umap", label = T,group.by = "celltype.3rd")
(p|p1) / (p2|p3)

# 保存结果
saveRDS(sc_data,"Seurat_celltype.rds")





#### 写在后面 ####
# 依我愚见，除非你是在构建精细的细胞图谱（Cell Atlas），
# 否则在绝大多数机制研究或临床样本对比中，
# 并不需要强求将每一个Cluster、每一个细胞类型都注释得清清楚楚。
# 这往往会陷入“为了注释而注释”的陷阱，投入产出比极低，并且客观上，无法将所有细胞注释清楚。

# 更好的科研策略是：
# 1. 先快速完成大类注释，确保主要谱系（免疫、上皮、基质等）无误；
# 2. 锁定你感兴趣的、或者在生物学上有显著差异的目标大类；
# 3. 提取（Subset）出来，重新跑降维聚类。
# 在去除了其他细胞干扰的环境下，针对性地进行精细注释。
# 简而言之：抓大放小，聚焦核心。

#### 亚群提取 ####
# 参数设置
# 要使用的细胞类型列名
celltype_col <- "celltype.3rd"
# 要提取的细胞大类
target_celltype <- "Excitatory"

# 提取Seurat对象
sc_data.sub <- sc_data[, sc_data@meta.data[[celltype_col]] %in% target_celltype]

# 保存结果
saveRDS(sc_data.sub,"Seurat_celltype_Excitatory.rds")

#### 快速演示亚型注释 ####
# 整体步骤：重新标准化-降维-Harmony后，再进行细胞类型注释

# 加载R包
library(harmony)

# 标准化
sc_data.sub = sc_data.sub %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()  

# PCA
sc_data.sub <- RunPCA(sc_data.sub)
pcs = 1:30

# Harmony
sc_data.sub <- RunHarmony(object = sc_data.sub,  
                          group.by.vars = "orig.ident", # 指定批次变量
                          reduction = "pca",            # 基于PCA结果校正  
                          dims = pcs,                   # 使用选定的PC数
                          assay.use="RNA",              # 选择特定的assay，默认
                          plot_convergence = TRUE,      # 显示收敛曲线
                          max.iter = 20,                # 最大迭代次数
                          reduction.save = "harmony"    # 保存结果的名称为"harmony"   
)

# UMAP和t-SNE
sc_data.sub <- FindNeighbors(sc_data.sub,
                             reduction = "harmony",  
                             dims = pcs) %>% 
               FindClusters(resolution = 0.1) %>% 
               RunUMAP(reduction = "harmony",  
                       dims = pcs) %>% 
               RunTSNE(dims = pcs, 
                       reduction = "harmony")

DimPlot(sc_data.sub, reduction = "umap", label = T)
DimPlot(sc_data.sub, reduction = "umap", label = F, group.by = "orig.ident")
DimPlot(sc_data.sub, reduction = "tsne", label = T)
DimPlot(sc_data.sub, reduction = "tsne", label = F, group.by = "orig.ident")

# 保存结果
saveRDS(sc_data.sub,"Seurat_celltype_Excitatory_harmony_30_0.1.rds")