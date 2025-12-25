# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面
# 基本原理是基于 G2/M 和 S 期的 markers 基因表达情况
# 利用CellCycleScoring函数，对各细胞进行打分，判断其所处阶段（phase）
# 本部分为可选步骤，并非所有单细胞分析中都要进行

# 这一步在读取数据且进行过QC后方可运行
# 运行完本部分后，接下来为PCA、降维聚类

#### 加载数据 ####
# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

# 设置工作目录
setwd("文件目录")

# 读取已经经过初步质控的Seurat对象
sc_data = readRDS("Seurat_qc.rds")

# Seurat 包中内置了一组和细胞周期相关的markers
str(cc.genes)      # Seurat内置的细胞周期markers
cc.genes$s.genes   # S期markers,43个
cc.genes$g2m.genes # G2M期markers,54个

# 获取G2M期相关基因，并与数据进行匹配
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(sc_data))

# 获取S期相关基因，并与数据进行匹配
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(sc_data))

#### 仅作拓展用 ####
# 但是仅针对人类这一物种
# 故如果分析其他物种，需要使用biomaRt等工具进行同源匹配。
# https://github.com/hbc/tinyatlas/tree/master/cell_cycle
# 包含：Danio_rerio、Drosophila_melanogaster、Homo_sapiens、Mus_musculus
# 参考资料：https://mp.weixin.qq.com/s/yanrXoGox9kpH9-FIGQSwA

#### 核心代码 ####
##### 1.Normalize+FindVariableFeatures+ScaleData #####
# 1.标准化Seurat对象
sc_data = sc_data %>%
  NormalizeData() %>%         # 对UMI计数数据进行标准化
  FindVariableFeatures() %>%  # 寻找高变基因用于下游分析
  ScaleData()                 # 对基因表达进行缩放（均值为0，方差为1）

# 2.细胞周期阶段评分
sc_data <- CellCycleScoring(sc_data,
                            g2m.features=g2m_genes,
                            s.features=s_genes)
# 3.回归细胞周期的影响
sc_data = ScaleData(sc_data, 
                    vars.to.regress = c("S.Score", "G2M.Score"))

# 保存去除细胞周期后的结果
saveRDS(sc_data,"Seurat_CellCycle.rds")

##### 2.SCTransform #####
# 这部分没有标准答案，存在两种观点：
# https://github.com/satijalab/seurat/issues/1679
# 1.先进行Normalize，随后CellCycleScoring，最后SCTransform
# 2.先进行SCTransform(RNA assay)，随后CellCycleScoring，最后SCTransform(SCT assay)
# 好在似乎对于结果的影响不大，但甚至是剔除之后在可视化仍然分散(打咩！)
# 因此我选择代码更简洁，运算速度更快的方案1

# 1.Normalize
# 经查阅，CellCycleScoring函数使用data数据进行运算，所以应该只进行NormalizeData即可
sc_data = NormalizeData(sc_data) # 对UMI计数数据进行标准化

# 2.细胞周期阶段评分
sc_data <- CellCycleScoring(sc_data,
                            g2m.features=g2m_genes,
                            s.features=s_genes)

# 3.SCTransform（标准化的同时，回归细胞周期）
sc_data <- SCTransform(
  sc_data,
  vars.to.regress = c("S.Score", "G2M.Score"),
  verbose = FALSE
)

# 保存去除细胞周期后的结果
saveRDS(sc_data,"Seurat_CellCycle.rds")

#### 可视化讲解 ####
##### 1.Normalize+FindVariableFeatures+ScaleData #####
# 标准化Seurat对象
sc_data = sc_data %>%
  NormalizeData() %>%         # 对UMI计数数据进行标准化
  FindVariableFeatures() %>%  # 寻找高变基因用于下游分析
  ScaleData()                 # 对基因表达进行缩放（均值为0，方差为1）

# 计算细胞周期阶段评分 
sc_data <- CellCycleScoring(sc_data,
                            g2m.features=g2m_genes,
                            s.features=s_genes)

# 原始数据可视化
# PCA
sc_data <- RunPCA(sc_data)
pcs = 1:40

# 基于PCA的可视化
p_rna_pca_before <- DimPlot(sc_data, reduction = "pca", group.by = "Phase") + 
  ggtitle("Before cell cycle regression")
# p_rna_pca_before

# 降维聚类
sc_data <- FindNeighbors(sc_data, 
                         reduction = "pca",  
                         dims = pcs) %>% 
           FindClusters(resolution = 0.5) # 根据情况选择分辨率

sc_data <- RunUMAP(sc_data, 
                   reduction = "pca",  
                   dims = pcs) %>% 
           RunTSNE(dims = pcs, 
                   reduction = "pca")

# 基于UMAP/t-SNE的可视化
p_rna_umap_before <- DimPlot(sc_data,
                             reduction = "umap",
                             group.by = "Phase") +
  ggtitle("UMAP - Before cell cycle regression")
p_rna_tsne_before <- DimPlot(sc_data,
                             reduction = "tsne",
                             group.by = "Phase") +
  ggtitle("t-SNE - Before cell cycle regression")

# 回归细胞周期的影响
sc_data = ScaleData(sc_data, 
                    vars.to.regress = c("S.Score", "G2M.Score"))

# 剔除细胞周期影响后可视化
# PCA
sc_data <- RunPCA(sc_data)
pcs = 1:40

# 基于PCA的可视化
p_rna_pca_after <- DimPlot(sc_data, reduction = "pca", group.by = "Phase") + 
  ggtitle("After cell cycle regression")

# 降维聚类
sc_data <- FindNeighbors(sc_data, 
                         reduction = "pca",  
                         dims = pcs) %>% 
  FindClusters(resolution = 0.5) # 根据情况选择分辨率

sc_data <- RunUMAP(sc_data, 
                   reduction = "pca",  
                   dims = pcs) %>% 
  RunTSNE(dims = pcs, 
          reduction = "pca")

# 基于UMAP/t-SNE的可视化
p_rna_umap_after <- DimPlot(sc_data,
                             reduction = "umap",
                             group.by = "Phase") +
  ggtitle("UMAP - After cell cycle regression")
p_rna_tsne_after <- DimPlot(sc_data,
                             reduction = "tsne",
                             group.by = "Phase") +
  ggtitle("t-SNE - After cell cycle regression")

# 可视化对比
p_rna_pca_before | p_rna_pca_after
p_rna_umap_before | p_rna_umap_after
p_rna_tsne_before | p_rna_tsne_after

# 保存去除细胞周期后的结果
saveRDS(sc_data,"Seurat_CellCycle.rds")

##### 2.SCTransform #####
# 标准化Seurat对象
sc_data = sc_data %>%
  NormalizeData() %>%         # 对UMI计数数据进行标准化
  FindVariableFeatures() %>%  # 寻找高变基因用于下游分析
  ScaleData()                 # 对基因表达进行缩放（均值为0，方差为1）

# 计算细胞周期阶段评分 
sc_data <- CellCycleScoring(sc_data,
                            g2m.features=g2m_genes,
                            s.features=s_genes)

# 原始数据可视化
# PCA
sc_data <- RunPCA(sc_data)
pcs = 1:40

# 基于PCA的可视化
p_rna_pca_before <- DimPlot(sc_data, reduction = "pca", group.by = "Phase") + 
  ggtitle("Before cell cycle regression")
# p_rna_pca_before

# 降维聚类
sc_data <- FindNeighbors(sc_data, 
                         reduction = "pca",  
                         dims = pcs) %>% 
  FindClusters(resolution = 0.5) # 根据情况选择分辨率

sc_data <- RunUMAP(sc_data, 
                   reduction = "pca",  
                   dims = pcs) %>% 
  RunTSNE(dims = pcs, 
          reduction = "pca")

# 基于UMAP/t-SNE的可视化
p_rna_umap_before <- DimPlot(sc_data,
                             reduction = "umap",
                             group.by = "Phase") +
  ggtitle("UMAP - Before cell cycle regression")
p_rna_tsne_before <- DimPlot(sc_data,
                             reduction = "tsne",
                             group.by = "Phase") +
  ggtitle("t-SNE - Before cell cycle regression")

# 回归细胞周期的影响
# 加速SCT运算的R包
# BiocManager::install("glmGamPoi")
sc_data <- SCTransform(
  sc_data,
  vars.to.regress = c("S.Score", "G2M.Score"),
  verbose = FALSE
)

# 剔除细胞周期影响后可视化
# PCA
sc_data <- RunPCA(sc_data)
pcs = 1:40

# 基于PCA的可视化
p_sct_pca_after <- DimPlot(sc_data, reduction = "pca", group.by = "Phase") + 
  ggtitle("After cell cycle regression")

# 降维聚类
sc_data <- FindNeighbors(sc_data, 
                         reduction = "pca",  
                         dims = pcs) %>% 
  FindClusters(resolution = 0.5) # 根据情况选择分辨率

sc_data <- RunUMAP(sc_data, 
                   reduction = "pca",  
                   dims = pcs) %>% 
  RunTSNE(dims = pcs, 
          reduction = "pca")

# 基于UMAP/t-SNE的可视化
p_sct_umap_after <- DimPlot(sc_data,
                            reduction = "umap",
                            group.by = "Phase") +
  ggtitle("UMAP - After cell cycle regression")
p_sct_tsne_after <- DimPlot(sc_data,
                            reduction = "tsne",
                            group.by = "Phase") +
  ggtitle("t-SNE - After cell cycle regression")

# 可视化对比
p_rna_pca_before | p_sct_pca_after
p_rna_umap_before | p_sct_umap_after
p_rna_tsne_before | p_sct_tsne_after

# 保存去除细胞周期后的结果
saveRDS(sc_data,"Seurat_CellCycle.rds")

#### 补充方案 ####
# 原理说明：
# 1.直接回归 S.Score + G2M.Score 会移除所有与细胞周期相关的信号
#   在发育/分化/干细胞等场景下，"是否在增殖" 本身可能是重要生物学信息
# 2.通过仅回归 S.Score - G2M.Score（CC.Difference）
#   保留：非增殖（G1） vs 增殖（S/G2M）之间的差异
#   去除：增殖细胞内部 S 期与 G2M 期的阶段性差异
# 3.该策略适用于细胞分化、发育轨迹等分析场景

# 官方举例：研究小鼠造血，干细胞处于静止状态而分化细胞处于增殖状态（反之亦然）
# 去除所有细胞周期效应也会模糊干细胞和祖细胞之间的界限

# 代码演示
# 标准化 Seurat 对象
sc_data = sc_data %>%
  NormalizeData() %>%         # 对 UMI 计数数据进行标准化
  FindVariableFeatures()      # 寻找高变基因用于下游分析

# 计算细胞周期阶段评分
# S.Score 和 G2M.Score 分别反映 S 期和 G2/M 期相关基因的表达水平
sc_data <- CellCycleScoring(sc_data,
                            g2m.features=g2m_genes,
                            s.features=s_genes)

# 构建细胞周期差值指标
# CC.Difference 主要刻画：增殖细胞中 S 期 vs G2M 期的差异
sc_data$CC.Difference <- sc_data$S.Score - sc_data$G2M.Score

# 回归 CC.Difference，仅移除增殖细胞内部的细胞周期阶段效应
sc_data <- ScaleData(
  sc_data,
  vars.to.regress = "CC.Difference"
)

# 可视化部分省略
# 保存结果
saveRDS(sc_data, "Seurat_CellCycle_CCdifference.rds")
