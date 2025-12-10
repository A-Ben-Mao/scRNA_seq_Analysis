# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

#### 写在前面 ####
# 初步质控后的计数矩阵仍然存在技术偏差，最主要的是测序深度差异
# 即每个细胞的总测序量不同
# 标准化的目的是消除这种差异，使细胞间的表达量具有可比性
# 但需要注意的是：标准化≠批次校正

# 加载R包
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(patchwork)
library(sctransform)

# 设置工作目录
setwd("文件目录")

# 读取已经经过初步质控的Seurat对象
sc_data = readRDS("Seurat_qc.rds")

#### Normalizedata, Findvariablefeatures和Scaledata ####
# %>%（管道操作符，forward-pipe operator）是最常用的一种操作符，
# 就是把左侧准备的数据或表达式，传递给右侧的函数调用或表达式进行运行，
# 可以连续操作就像一个链条一样。

# 简写版本
sc_data = sc_data %>%
  NormalizeData() %>%         # 对UMI计数数据进行标准化
  FindVariableFeatures() %>%  # 寻找高变基因用于下游分析
  ScaleData()                 # 对基因表达进行缩放（均值为0，方差为1）
                              # 使得高表达基因不会占主导地位

# 可调整参数版本
sc_data <- NormalizeData(sc_data, 
                         normalization.method = "LogNormalize", 
                         scale.factor = 10000)
sc_data <- FindVariableFeatures(sc_data, 
                                selection.method = "vst", 
                                nfeatures = 2000)
sc_data <- ScaleData(sc_data) # 仅针对高变基因进行放缩，Seurat官方推荐

# 如果有特殊需求，可根据全基因进行放缩
# 同时可以回归掉一些协变量所带来的影响
# sc_data <- ScaleData(sc_data,
#                      vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
#                      features = rownames(sc_data)) # 全基因放缩

# 完整版本
?NormalizeData
?FindVariableFeatures
?ScaleData

# 默认的矩阵，可以根据情况切换
DefaultAssay(sc_data)
DefaultAssay(sc_data) = "RNA"
DefaultAssay(sc_data)

# 高变基因的查看与可视化
plot1 <- VariableFeaturePlot(sc_data)
rna_top10 <- head(VariableFeatures(sc_data), 10)
plot2 <- LabelPoints(plot = plot1, points = rna_top10)
plot1 | plot2  # 并排显示

#### 保存完成标准化的文件 ####
saveRDS(sc_data,"Seurat_Normalize.rds")

#### SCTransform ####
# 相当于替代了上述的三个函数NormalizeData，FindVariable，ScaleData
# Seurat官方似乎更加推荐此方法
# 对测序深度的校正效果要好于log标准化(10万以内的细胞都建议使用SCT)
# 寻找3000个高变基因
# 相对的，该R包对配置（特别是内存）要求更高，且运行时间更长

# 官方声明可以使用这个R包加速运行
# 如果已安装，则默认使用该包
# install.packages('BiocManager')
# BiocManager::install('glmGamPoi')

# 结果存储在sc_data@assays$SCT中
?SCTransform
sc_data <- SCTransform(sc_data)

# 同样的此方法也可以回归掉一些协变量所带来的影响
# sc_data <- SCTransform(sc_data, vars.to.regress = "percent.mt")
# sc_data <- SCTransform(sc_data, vars.to.regress = c("S.Score", "G2M.Score"))

# 默认的矩阵，可以根据情况切换
DefaultAssay(sc_data)
DefaultAssay(sc_data) = "SCT"
DefaultAssay(sc_data)

# 高变基因的查看与可视化
plot3 <- VariableFeaturePlot(sc_data)
sct_top10 <- head(VariableFeatures(sc_data), 10)
plot4 <- LabelPoints(plot = plot3, points = sct_top10)
plot3 | plot4  # 并排显示

#### 保存完成标准化的文件 ####
saveRDS(sc_data,"Seurat_Normalize_SCTchangeRNA.rds")

# 需要注意的是
# 查阅相关资料得知，两种标准化方法的结果是存放在不同的层级中，理应相互独立
# 但后续实际操作的分析中会发现似乎后进行的标准化方法会“抵消”另一部分的结果
# 比如：若先进行RNA标准化，后进行SCT标准化时，高变基因的信息似乎会丢失
# 不知是否和R包的版本相关
# 若确实想用一个Seurat对象保存两种标准化，建议先SCT后切换assay进行Normalize
# 否则推荐选择一种标准化方法进行全流程

