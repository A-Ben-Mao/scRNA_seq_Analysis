# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 这一步需要建立在已经整合读取多个样本的数据之后

# 设置工作目录
setwd("文件目录")

# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

# 加载已经整合读取的Seurat对象
sc_data <- readRDS("Seurat_combined.rds")

# 查看样本来源
table(sc_data@meta.data$orig.ident)

#### 质量控制 ####
# 线粒体基因比例
# 高比例可能表示细胞死亡或损伤
# 线粒体基因的命名方式，通常以“MT-”开头（人类）或“mt-”开头（小鼠）
sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")

# 红细胞比例（可选）
# 高比例可能表示红细胞污染
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(sc_data)) # 匹配实际存在的基因
sc_data[["percent.HB"]]<-PercentageFeatureSet(sc_data, features=HB.genes) 

# 质量评估可视化
# 查看相关性
FeatureScatter(sc_data, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(sc_data, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

# 绘制质控前小提琴图
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(sc_data, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
violin

# 保存图像
# ggsave("data_violin_before_qc.pdf", plot = violin, width = 14, height = 8) 
dim(sc_data)

# 查看数据分布特征
# Feature
quantile(sc_data$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(sc_data$nFeature_RNA, seq(0.9, 1, 0.01))
# Count
quantile(sc_data$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(sc_data$nCount_RNA, seq(0.9, 1, 0.01))
# 线粒体比例
quantile(sc_data$percent.mt, seq(0.9, 1, 0.01))
# 红细胞比例
quantile(sc_data$percent.HB, seq(0.9, 1, 0.01))

# 设置质控标准
# Feature
minGene=200
maxGene=10000 # 或许上限设置较低可免除双细胞剔除
              # 如果有能力，还是建议根据尝试双细胞之后看那种结果更合适
              # 当然也可以不进行双细胞的剔除
# Count
minUMI=200
# 线粒体比例
pctMT=10
# 血细胞比例
pctHB=1

# 不同文献的阈值都有所不同，大家根据自己数据集所在原文献的标准进行质控即可

# 应用质控过滤
sc_data <- subset(sc_data, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)

# 绘制质控后小提琴图（验证过滤效果）
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(sc_data, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
violin

# 保存图像
# ggsave("data_violin_after_qc.pdf", plot = violin, width = 14, height = 8) 
dim(sc_data)

# 保存质控后Seurat对象
saveRDS(sc_data,"Seurat_qc.rds")

