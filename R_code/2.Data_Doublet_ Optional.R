# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面
# DoubletFinder官方流程推荐单样本处理
# 每个单样本Seurat对象至少需要进行到PCA，本代码进行了细胞聚类以可视化呈现

# 同时，因为该方法对同型双联体（即由转录相似的细胞状态形成的双联体）不敏感
# 故而建议利用文献支持的细胞类型注释来模拟数据中存在的同型双联体的比例

# 本示例使用GSE302285，为10X数据
# 其他数据类型同理，可借用AI修改代码
# 整体原则：在读取数据后的merge和joinlayer之前，提取list中的对象进行去除双细胞

# 完成这一步后再进行标准化，以及后续操作

# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)
library(ggplot2)
library(patchwork)

##### matrix.mtx、genes.tsv和barcodes.tsv 的读取####
# 需要读取的文件目录
setwd("文件目录")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files()
samples

# 创建一个空的列表
seurat_list <- list()

#读取数据并创建Seurat对象
#删除，小于200个基因表达的细胞，小于3个细胞表达的基因
for (sample in samples) {
  
  #读取10x数据
  seurat_data <- Read10X(data.dir = sample)
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = seurat_data, # 原始表达矩阵
    project = sample,     # 为对象设置项目名，metadata中可区分来源
    min.features = 200,   # 细胞阈值
    min.cells = 3         # 基因阈值
  )
  
  #添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

# 删除中间数据，防止内存溢出
rm(seurat_data)
rm(seurat_obj)

#### 提取样本 ####
# 由于此操作易致内存溢出，故推荐逐个进行，而非循环处理
# 同时不同细胞量对应的双细胞比例不同（见后）
case <- seurat_list[[1]] # 修改不同样本的序号

#### 单样本质控 ####
# 线粒体基因比例，“MT-”开头（人类）或“mt-”开头（小鼠）
case[["percent.mt"]] <- PercentageFeatureSet(case, pattern = "^MT-")

# 红细胞比例（可选）
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(case)) # 匹配实际存在的基因
case[["percent.HB"]]<-PercentageFeatureSet(case, features=HB.genes) 

# 质量评估可视化
# 查看相关性
FeatureScatter(case, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(case, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

# 绘制质控前小提琴图
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(case, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
violin

# 查看数据分布特征
# Feature
quantile(case$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(case$nFeature_RNA, seq(0.9, 1, 0.01))
# Count
quantile(case$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(case$nCount_RNA, seq(0.9, 1, 0.01))
# 线粒体比例
quantile(case$percent.mt, seq(0.9, 1, 0.01))
# 红细胞比例
quantile(case$percent.HB, seq(0.9, 1, 0.01))

# 设置质控标准
# Feature
minGene=200
maxGene=10000 
# Count
minUMI=200
# 线粒体比例
pctMT=10
# 血细胞比例
pctHB=1

# 应用质控过滤
case <- subset(case, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                    nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)

# 绘制质控后小提琴图（验证过滤效果）
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(case, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
violin

#### DoubletFinder双细胞识别 ####
# DoubletFinder可分为4个步骤:
# (1)从已有的scRNA-seq数据中生成人工doublet;
# (2)对合并的真实人工数据进行预处理;
# (3)进行主成分分析，利用PC距离矩阵求出每个单元的人工k个最近邻(pANN)的比例;
# (4)根据期望的doublet数量排序和阈值pANN值。

# 标准化(也可选择SCTransform)
case = case %>%
  NormalizeData() %>%         # 对UMI计数数据进行标准化
  FindVariableFeatures() %>%  # 寻找高变基因用于下游分析
  ScaleData()                 # 对基因表达进行缩放（均值为0，方差为1）

# PCA降维
case <- RunPCA(case)

# 个性化选择参数(可选)，直接默认1:30也可
# 画出前 50 个 PC 的标准差或方差贡献
# ElbowPlot(case, ndims = 50)
# pct <- case [["pca"]]@stdev / sum( case [["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# cumu

# PCA、聚类、非线形降维
case = case %>% 
  RunUMAP(dims = 1:30) %>%        # 非线性降维可视化
  RunTSNE(dims = 1:30) %>%        # 非线性降维可视化
  FindNeighbors(dims = 1:30) %>%  # 构建细胞间的邻域图
  FindClusters(resolution = 0.1)  # 基于图论的细胞聚类（分辨率0.1得到较少的聚类）

# 扫描不同pK（邻域大小）参数(运行时间较长)
sweep.res.list <- paramSweep(case, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# 选择最佳pK
best_pk <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 估算同型双细胞比例（难以检测）
annotations <- case$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# 估算异型双细胞（相对容易检测）
# 最终测量不同细胞数，其双细胞比例不同，可参考：
# https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
dim(case) # 数据维度(基因数 细胞数)
nExp_poi <- round(0.075*nrow(case@meta.data)) # 默认比例为0.075

# 也有前辈进行动态计算预期双细胞率 (线性估算: 每1000个细胞增加约0.8%的双细胞率)
# n_cells <- nrow(case@meta.data)
# expected_rate <- (n_cells / 1000) * 0.008 
# nExp_poi <- round(expected_rate * n_cells)

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))  # 调整预测数(排除同源双细胞干扰)

# 运行DoubletFinder
case <- doubletFinder(case, 
                      PCs = 1:30, 
                      pN = 0.25,           # 人工双细胞比例(默认0.25)
                      pK = best_pk,        # 之前获得的最优邻居比例参数
                      nExp = nExp_poi.adj, # nExp: 预期双细胞数量
                      reuse.pANN = NULL,   # 首次运行计算双细胞概率结果
                      sct = FALSE)

# 将列名改为"Doublet_score"和"Singlet_Doublet"
colnames(case@meta.data)
colnames(case@meta.data)[length(colnames(case@meta.data))-1] <- "Doublet_score"
colnames(case@meta.data)[length(colnames(case@meta.data))] <- "Singlet_Doublet"
colnames(case@meta.data)

# 查看DoubletFinder分析结果
head(case@meta.data[, c("Doublet_score", "Singlet_Doublet")])

# 可视化
# 绘制DoubletFinder分类的umap/tsne图
DimPlot(case, reduction = "umap", group.by = "Singlet_Doublet")
DimPlot(case, reduction = "tsne", group.by = "Singlet_Doublet")

# 绘制双细胞分类的小提琴图
VlnPlot(case, group.by = "Singlet_Doublet", 
        features = c("nCount_RNA", "nFeature_RNA"), 
        pt.size = 0.1, 
        ncol = 2)

# 过滤非单细胞数据
case <- subset(case, Singlet_Doublet == "Singlet")
cat("Cells after filtering singlets:", ncol(case), "\n")

# 清理冗余数据(似乎不够彻底，但后续分析会覆盖)
if(length(names(case@reductions)) > 0){
  for(rn in names(case@reductions)) case[[rn]] <- NULL
}

# 将处理好的Seurat对象返回list
seurat_list[[1]] <- case # 修改不同样本的序号

#### 合并处理好的数据 ####
# 合并多个 Seurat 对象（跨样本合并）
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

# 将Layers融合
sc_data = JoinLayers(seurat_combined)

# 保存Seurat对象

saveRDS(sc_data, file = "Seurat_combined.rds")
