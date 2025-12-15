# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面
# 本步骤需要在进行了质控和标准化处理后运行

# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(clustree)

# 设置工作目录
setwd("文件目录")

# 读取已经经过标准化的Seurat对象
sc_data = readRDS("Seurat_Normalize.rds")

# 设置随机数(涉及UMAP/TSNE降维)
set.seed(1234)

# 查看对象中所有的assay
print(sc_data@assays)

# 查看当前选择的assay，选择特定的assay
DefaultAssay(sc_data)
DefaultAssay(sc_data) = "SCT"
DefaultAssay(sc_data) = "RNA"
DefaultAssay(sc_data)

# 返回当前细胞的身份分组
# 后续作图可根据不同分组进行探索
# 如可以探索单双细胞、不同细胞周期、不同Cluster
table(Idents(sc_data))
# Idents(sc_data) = "orig.ident"

####  PCA Principal Component Analysis ####
# 结果存进sc_data@reductions$pca
sc_data <- RunPCA(sc_data)

# 可视化 PCA 的前两个维度（PC1 vs PC2）的细胞分布
# 通常用于检查样本/群体是否分离或是否存在批次效应
DimPlot(object = sc_data, reduction = "pca")

# 绘制每个PCA成分的相关基因
# 显示每个 PC 中对该主成分贡献最大的基因
VizDimLoadings(object = sc_data, dims = 1:6, reduction = "pca",nfeatures = 20)

# 绘制与每个 PC 相关的细胞/基因热图
DimHeatmap(sc_data, dims = 1:6, cells = 500, balanced = TRUE)

# 画出前 50 个 PC 的标准差或方差贡献
ElbowPlot(sc_data, ndims = 50)

# 确定与每个 PC 的百分比并计算累计百分比
pct <- sc_data [["pca"]]@stdev / sum( sc_data [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
cumu 

# 选择合适的PC
# 通常选择20-40个PC
# 一般选择累积贡献达到 80%/90% 时的 PC 数
pcs = 1:40 

#### 细胞聚类 Cell Clustering ####
# 目的：选取合适的分辨率
# 基于指定的 PC 计算邻接图（kNN graph）
sc_data <- FindNeighbors(sc_data,  dims = pcs) 

# 从0.1-2的resolution结果均运行一遍
seq = seq(0.1,2,by=0.1)
for (res in seq){
  sc_data = FindClusters(sc_data, resolution = res)
}

# 绘制树状图
p1 = clustree(sc_data,prefix = "RNA_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
p

#### 非线性降维 (UMAP/tSNE) ####
# 降维聚类
sc_data <- FindNeighbors(sc_data, 
                         reduction = "pca",  
                         dims = pcs) %>% 
           FindClusters(resolution = 0.5) # 根据情况选择分辨率

sc_data <- RunUMAP(sc_data, 
                   reduction = "pca",  
                   dims = pcs) %>% 
           RunTSNE(dims = pcs, 
                   # check_duplicates = FALSE, # 若因重复而报错，则添加
                   reduction = "pca")

# 可视化
DimPlot(sc_data, reduction = "umap", label = T)

DimPlot(sc_data, reduction = "umap", label = F, group.by = "orig.ident")

DimPlot(sc_data, reduction = "tsne", label = T)

DimPlot(sc_data, reduction = "tsne", label = F, group.by = "orig.ident")

dev.off()

# 保存结果
saveRDS(sc_data,"Seurat_DimensionalityReduction.rds")

