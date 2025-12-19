# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面
# SingleR自带7个参考数据集：
# https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html
# 5个是人类数据
# BlueprintEncodeData Blueprint (Martens and Stunnenberg 2013) and Encode (The ENCODE Project Consortium 2012) （人）
# DatabaseImmuneCellExpressionData The Database for Immune Cell Expression(/eQTLs/Epigenomics)(Schmiedel et al. 2018)（人）
# HumanPrimaryCellAtlasData the Human Primary Cell Atlas (Mabbott et al. 2013)（人）
# MonacoImmuneData, Monaco Immune Cell Data - GSE107011 (Monaco et al. 2019)（人）
# NovershternHematopoieticData Novershtern Hematopoietic Cell Data - GSE24759（人）
# 2个是小鼠的数据
# ImmGenData the murine ImmGen (Heng et al. 2008) （鼠）
# MouseRNAseqData a collection of mouse data sets downloaded from GEO (Benayoun et al. 2019).鼠）

# SingleR将这几个数据库单独放在了celldex  R包里
# 每次调用需要下载，故而推荐将参考数据集下载到本地存储，每次直接调用即可
# 相关网盘下载链接已放入“reference_data”文件夹中，若失效可根据下列代码自行下载

# BiocManager::install("celldex") # 如果需要自行下载参考数据集可下载
# library(celldex)
# ref <- HumanPrimaryCellAtlasData() # 14M
# ref <- BlueprintEncodeData() # 292M
# ref <- MouseRNAseqData() # 42M
# ref <- ImmGenData() # 56M
# ref <- DatabaseImmuneCellExpressionData() # 24M
# ref <- NovershternHematopoieticData() # 16M
# ref <- MonacoImmuneData() #8.6M

# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
# BiocManager::install("SingleR")
library(SingleR)

# 设置工作目录
setwd("文件目录")

# 读取已经经过降维的Seurat对象
sc_data = readRDS("Seurat_DimensionalityReduction_harmony_40_0.5.rds")

# 加载参考数据库
load("E:\\scRNA-seq Analysis\\SingleR_ref\\HumanPrimaryCellAtlasData.rda")

# 设置参考参考集名称
ref_data <- ref_hpca

# 查看参考数据集中包含的细胞类型
# 主要细胞类型（如T细胞、B细胞等）
table(ref_data@colData@listData[["label.main"]])
# 更精细的亚型分类
table(ref_data@colData@listData[["label.fine"]])

# 获取表达矩阵
testdata = GetAssayData(object = sc_data@assays$RNA, layer = "data")
# testdata = GetAssayData(object = sc_data@assays$SCT, layer = "data")

# 获取clusters信息
clusters <- sc_data@meta.data$seurat_clusters

# 运行SingleR
cellpred <- SingleR(test = testdata,
                    ref = ref_data,
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    labels = ref_data@colData@listData[["label.main"]], # 可修改为label.fine
                    assay.type.ref = "logcounts")

# 将聚类水平的注释结果映射到每个细胞
celltype = data.frame(ClusterID=rownames(cellpred),
                      celltype=cellpred$labels,
                      stringsAsFactors = F)

# 在metadata中创建新的SingleR列存储注释结果
sc_data@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
  sc_data@meta.data[which(sc_data$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
}

# 注释结果可视化
# 显示每个聚类与参考细胞类型的匹配得分
p = plotScoreHeatmap(cellpred)

# 显示t-SNE和UMAP结果
p1 <- DimPlot(sc_data, group.by = "SingleR", label = T,reduction = "tsne")
p2 <- DimPlot(sc_data, group.by = "SingleR", label = T,reduction = "umap")
p3 <- p1 | p2
p3

# 保存结果文件
saveRDS(sc_data,"Seurat_SingleR.rds")
