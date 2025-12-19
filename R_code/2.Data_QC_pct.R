# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 这一步需要建立在已经整合读取多个样本的数据之后
# 这是百分比版本的QC代码，可根据数据特征选择该方法进行质控

# Q:为什么采用“分样本百分比”而非“全局固定值”过滤？
# 1.异质性：不同样本的测序深度和细胞类型复杂度差异可能很大。
# 2.避免误删：若使用全局固定上限（如 nFeature < 5000），测序深度较浅的样本可能全是正常细胞，而深度极高的样本中可能混入了大量双细胞（Doublet）却未被剔除。
# 3.相对公平：使用百分位数（如 95%）能根据每个样本自身的分布特征，自适应地切除该样本中的离群值。

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

#### 设置质控标准 ####
# 1. 下限 (固定值，去除空液滴)
minGene = 200
minUMI  = 200  

# 2. 上限(百分比，分开设置，自适应去除双细胞)
# 设置 0.95 代表去除该样本中最高的 5%
high_quantile_gene = 0.95  # 针对 Gene (nFeature_RNA) 的百分比
high_quantile_umi  = 0.95  # 针对 UMI (nCount_RNA) 的百分比

# 3. 比例阈值(固定值，去除死细胞/红细胞)
pctMT = 10
pctHB = 1

#### 应用质控过滤 ####
# 获取所有样本ID
sample_ids <- unique(sc_data$orig.ident)
cells_to_keep <- c() # 用于存放通过质控的细胞名

# 循环每个样本，单独计算该样本的上限
for(sid in sample_ids){
  # 提取当前样本的数据
  current_meta <- sc_data@meta.data[sc_data$orig.ident == sid, ]
  
  # 计算该样本特定的上限阈值
  maxGene_sample <- quantile(current_meta$nFeature_RNA, probs = high_quantile_gene)
  maxUMI_sample  <- quantile(current_meta$nCount_RNA,   probs = high_quantile_umi)
  
  # 筛选符合条件的细胞 ID
  good_cells <- rownames(current_meta)[
    current_meta$nFeature_RNA > minGene & 
      current_meta$nFeature_RNA < maxGene_sample &
      current_meta$nCount_RNA   > minUMI & 
      current_meta$nCount_RNA   < maxUMI_sample &
      current_meta$percent.mt   < pctMT & 
      current_meta$percent.HB   < pctHB
  ]
  
  # 将好细胞加入名单
  cells_to_keep <- c(cells_to_keep, good_cells)
  
  # 打印每个样本算出来的阈值
  print(paste0("样本 ", sid, 
               " | Gene上限(Top ", (1-high_quantile_gene)*100, "%): ", round(maxGene_sample), 
               " | UMI上限(Top ",  (1-high_quantile_umi)*100, "%): ", round(maxUMI_sample)))
}

# 将过滤后Seurat对象赋予sc_data
sc_data <- subset(sc_data, cells = cells_to_keep)

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
