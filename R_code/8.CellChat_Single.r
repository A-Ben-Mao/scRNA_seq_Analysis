# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面
# CellChat官方网站：https://github.com/jinworks/CellChat

# 这一步需要用到：细胞的基因表达数据，以及细胞标签
# 因此，Seurat对象需要至少运行到完成细胞注释

# 加载R包
library(Seurat)
library(CellChat)
library(patchwork)

# 设置工作目录
setwd("文件目录")

# 读取已经经过细胞注释的Seurat对象
sc_data <- readRDS("E:\\scRNA-seq Analysis\\data\\Seurat_celltype.rds")

# 查看当前 Idents，以及有哪些细胞类型
head(Idents(sc_data))
table(Idents(sc_data))

# 设置细胞类型信息
Idents(sc_data) <- "celltype.3rd" # 输入细胞注释列的列名，注意修改

# 由于原始的数据较大，运算时间较长，现提取部分数据进行演示
dim(sc_data)

set.seed(123)
library(dplyr)
meta <- sc_data@meta.data
meta$cell_id <- rownames(meta)
sampled_cells <- meta %>%
  group_by(celltype.3rd) %>%
  slice_sample(n = 500) %>% # 每种类型随机抽500个，不够的就全取
  pull(cell_id)

# 生成子集
sc_data <- subset(sc_data, cells = sampled_cells)
dim(sc_data)
table(Idents(sc_data))

# 创建并切换到新工作目录
output_dir <- "CellChat"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)  # 切换工作目录到目标文件夹

#### 数据准备 ####
##### 提取 CellChat 所需的两个核心对象 #####
# 1. 提取表达矩阵
data.input <- sc_data[["RNA"]]$data

# 2. 提取细胞分组信息
labels <- Idents(sc_data)
meta <- data.frame(labels = labels, row.names = names(labels))

# （可选）若只分析某一条件 / 样本
# 提取特定列的特定行
# 比如这里是“group“列的所有“疾病”行
# cells.use <- rownames(sc_data@meta.data)[sc_data$group == "disease"]
# data.input <- data.input[, cells.use]
# meta <- meta[cells.use, ]

##### 创建 CellChat 对象 #####
# group.by 决定：CellChat 按什么细胞分组来推断通讯
cellchat <- createCellChat(
  object   = sc_data,
  group.by = "ident",   # 使用 Idents(sc_data)，其实就是不同的细胞类型
  assay    = "RNA"      # 使用 RNA assay 的 log-normalized data
)

##### 设置配体-受体相互作用数据库 #####
# CellChatDB 是什么
# 它是一个“先验知识库”：人工整理、文献支持
# 定义了：哪个 ligand，可以和哪个 receptor，以哪种方式进行细胞通讯

# 选择对应的CellChatDB
CellChatDB <- CellChatDB.human # 人类
# CellChatDB <- CellChatDB.mouse # 小鼠

# 了解数据库结构
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
# 包括“分泌信号”、“细胞外基质-受体”、“细胞间接触”和“非蛋白信号”

# 设定选用的库
# 下面是官方示例的代码，根据情况选择
# 1. 只使用 Secreted Signaling
CellChatDB.use <- subsetDB(
  CellChatDB,
  search = "Secreted Signaling",
  key = "annotation"
)

# 2. 只使用来自V1的 Secreted Signaling
# CellChatDB.use <- subsetDB(
#   CellChatDB, 
#   search = list(c("Secreted Signaling"), 
#                 c("CellChatDB v1")), 
#   key = c("annotation", "version"))

# 3. 使用所有除开“非蛋白信号”的信号通路
# CellChatDB.use <- subsetDB(CellChatDB)

# 4. 使用所有的信号通路
# 由于包含“非蛋白信号”，官方并不是很推荐
# CellChatDB.use <- CellChatDB

# 将数据库写入 CellChat 对象
cellchat@DB <- CellChatDB.use

##### 预处理表达数据 #####
# 根据电脑配置选择并行运算核心数，以加速运行
future::plan("multisession", workers = 6)

# 1. 剔除无关基因
cellchat <- subsetData(cellchat)

# 2. 鉴别过表达基因
# 计算每个细胞群中哪些基因是显著高表达的（与所有其他组相比）
# 结果存储：cellchat@var.features
cellchat <- identifyOverExpressedGenes(cellchat)

# 3. 鉴别过表达相互作用
# 基于过表达基因和选择的数据库，匹配 配体-受体对
# 【注意】这里可能会由于超过默认的内存分配（500MiB）而报错
# 【方案】将限制设置为更大内存
# 【说明】仅对当前的 R Session 有效，重启RStudio后需要重新运行这段代码
options(future.globals.maxSize = 5 * 1024 * 1024 * 1024)
cellchat <- identifyOverExpressedInteractions(cellchat)

# （补充）将基因表达数据投影到“蛋白质-蛋白质相互作用（PPI）网络”
# 如果你的数据测序深度很低，这个功能可以填补配体或受体亚基的“假零值”
# 默认不进行，若初次运行结果较少，或者没有满意的相互作用，可以补充运行
# cellchat <- projectData(cellchat, PPI.human)

#### 推断细胞间通讯网络 ####
##### 计算通信概率并推断蜂窝通信网络 #####
# 1. 计算通讯概率（计算时间非常长）
# 默认方法为Trimean，基因在一个细胞群里的表达比例需要大于 25%才会被纳入分析
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# 注意：如果你之前做了 PPI 投影，这里 raw.use 要设为 FALSE

# 保存CellChat对象
saveRDS(cellchat, file = "cellchat_compute.rds")

# 加载已经计算的对象
# cellchat <- readRDS("cellchat_compute.rds")

# （拓展）使用 10% 截断均值，降低门槛
# 默认不进行，若初次运行结果较少，或者没有满意的相互作用，可以补充运行
# cellchat <- computeCommunProb(cellchat, 
#                               type = "truncatedMean", 
#                               trim = 0.1, # 可以尝试 0.1 或 0.05
#                               raw.use = TRUE)

# 2. 过滤低质量通讯，排除微小cluster的影响
cellchat <- filterCommunication(cellchat, min.cells = 10)

##### 提取结果 #####
# 下列代码为一个函数的不同用法，旨在提取不同信息
# 1. 提取所有推断出的配体-受体（Ligand-Receptor）层面的通讯结果（默认用法）
df.net <- subsetCommunication(cellchat)

# 2. 仅提取通路相关信息，而不关注LR信息
df.netP <- subsetCommunication(cellchat, slot.name = "netP")

# 3. 提取制定细胞群的结果
# levels(cellchat@idents)
df.net <- subsetCommunication(
  cellchat,
  sources.use = c(1,2), # 从细胞群1、2发出
  targets.use = c(4,5)  # 由细胞群4、5接收
)

df.net <- subsetCommunication(
  cellchat,
  sources.use = c("Oligodendrocytes", "Astrocytes"),
  targets.use = c("Inhibitory", "Pericytes")
)

# 4. 提取指定信号通路的结果
df.net <- subsetCommunication(
  cellchat,
  signaling = c("WNT", "TGFb")
)

##### 推断信号通路水平的细胞间通讯 #####
# 由于一些不同的LR对可能都是同一条信号通路，比如CCL
# 因此这一步是将同一通路的所有LR对的概率进行“汇总”
# 计算该通路的整体通讯概率
cellchat <- computeCommunProbPathway(cellchat)

##### 聚合的细胞间通信网络 #####
# 由于原始结果的数量级较大，且复杂，不易解读
# 因此通过这一步将大量的信息聚合成两个指标
# 数量：两个细胞群之间有多少对显著的配体-受体相互作用（cellchat@net$count）
# 强度：两个细胞群之间通讯的总概率之和是多少（cellchat@net$weight）
cellchat <- aggregateNet(cellchat)

##### 细胞通讯初步可视化 #####
# 1. 总体概览
# 计算每个细胞群的大小（细胞数量），用于决定圈图上节点的粗细
groupSize <- as.numeric(table(cellchat@idents))

# 设置画布为 1行2列
par(mfrow = c(1,2), xpd=TRUE)

# 左图：展示交互数量 (Number of interactions)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

# 右图：展示交互强度 (Interaction weights/strength)
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

# 2. 拆分图
mat <- cellchat@net$weight # 取出权重矩阵
par(mfrow = c(3,4), xpd=TRUE) # 设置画布布局（根据你的细胞群数量调整行列，这里是3行4列）

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

#### 细胞间通讯网络的可视化 ####
##### 信号通路层面的可视化 #####
# 1. 选定要展示的通路
pathways.show <- c("TGFb") 

# 2. 层次图 (Hierarchy Plot)
# vertex.receiver = seq(1,4) 意味着第 1,2,3,4 号细胞群作为“主要接收者”显示在左侧
vertex.receiver = seq(1,4) 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# 3. 圈图 (Circle Plot)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# 4. 和弦图 (Chord Diagram)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# 5. 热图 (Heatmap)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# （补充）分组和弦图 (netVisual_chord_cell)
# 如果你的细胞群分得很细（比如 Fibroblast_1, Fibroblast_2...），图会很乱
# 则可以将它们归类为大类（如 FIB, DC, TC）再画图
# 这里的数字“4”代表的是在 cellchat 对象的聚类顺序中，连续出现了多少个属于该大类的亚群
# levels(cellchat@idents)
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) 
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = ...)

##### 单个配体-受体对的可视化 #####
# 贡献度柱状图
netAnalysis_contribution(cellchat, signaling = pathways.show)

# 提取并绘制单个 L-R 对
# 1. 提取 CXCL 通路下所有显著的 L-R 对
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)

# 2. 选取第一对（例如 CXCL12-CXCR4）
LR.show <- pairLR.CXCL[1,] 

# 3. 针对这一对画图
# 层次图
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# 圈图
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

##### 自动循环绘图 #####
# 在初次分析时，若是手动逐个绘制太过耗费时间
# 因此采用循环代码自动逐个绘制
# 1. 获取所有显著通路的名称
pathways.show.all <- cellchat@netP$pathways

# 2. 设置接收者索引（为了层次图）
vertex.receiver = seq(1,4)

# 3. 循环遍历每一个通路
for (i in 1:length(pathways.show.all)) {
  
  # A. 画网络图 (这里使用了通用函数 netVisual，它会根据 layout 自动调用底层函数)
  # 注意：这里只画了层次图，你可以加代码画圈图
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  
  # B. 计算 L-R 贡献度并保存 PDF
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

#### 个性化的多配体-受体对或信号通路层面可视化 ####
##### 气泡图 #####
# 1. 指定细胞群之间的所有通讯
netVisual_bubble(cellchat, 
                 sources.use = 4, 
                 targets.use = c(5:11), 
                 remove.isolate = FALSE)

# 2. 特定细胞群之间的指定信号通路
netVisual_bubble(cellchat, 
                 sources.use = 4, 
                 targets.use = c(5:11), 
                 signaling = c("TGFb","WNT"), 
                 remove.isolate = FALSE)

# 3. 特定细胞群之间的指定LR对
# 提取感兴趣的LR对
pairLR.use <- extractEnrichedLR(cellchat, 
                                signaling = c("TGFb","WNT"))
# 图片绘制
netVisual_bubble(cellchat, 
                 sources.use = c(3,4), 
                 targets.use = c(5:8), 
                 pairLR.use = pairLR.use, 
                 remove.isolate = TRUE)

##### 和弦图 #####
# 1. 查看特定发送者 (Sender View)
netVisual_chord_gene(cellchat, 
                     sources.use = 4,  # 指定发送者
                     targets.use = c(5:11), 
                     lab.cex = 0.5,
                     legend.pos.y = 30)

# 2. 查看特定接收者 (Receiver View)
netVisual_chord_gene(cellchat, 
                     sources.use = c(1,2,3,4), 
                     targets.use = 8,  # 指定接收者
                     legend.pos.x = 15)

# 3. 查看通路层面的流向
netVisual_chord_gene(cellchat, 
                     sources.use = c(1,2,3,4), 
                     targets.use = c(5:11), 
                     signaling = c("TGFb","WNT"), # 指定通路
                     legend.pos.x = 8)

netVisual_chord_gene(cellchat, 
                     sources.use = c(1,2,3,4), 
                     targets.use = c(5:11), 
                     slot.name = "netP",  # 全部通路
                     legend.pos.x = 10)

##### 基因表达分布图（小提琴/点图） #####
# 默认用法 (enriched.only = TRUE)
plotGeneExpression(cellchat, 
                   signaling = "TGFb", 
                   enriched.only = TRUE, # 只画那些被 CellChat 识别为有显著通讯作用的配体和受体
                   type = "violin") # 可切换为dot、bar

#### 细胞间通讯网络的系统分析 ####
##### 识别细胞的信号角色 #####
# 1. 热图
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat,
                                  signaling = pathways.show,
                                  width = 8,
                                  height = 2.5,
                                  font.size = 10)

# 2. 散点图
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb", "WNT"))
gg1 + gg2

# 3. 信号贡献度热图
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

##### 识别并可视化分泌细胞的外向通讯模式 #####
# 加载R包
library(NMF)
library(ggalluvial)

# 确定模式数量
selectK(cellchat, pattern = "outgoing")

# 寻找 Cophenetic 和 Silhouette 指标突然下降的前一个点
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# 桑基图与点图可视化
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

##### 识别并可视化分泌细胞的传入通讯模式 #####
# 确定模式数量
selectK(cellchat, pattern = "incoming")

# 寻找参数
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# 桑基图与点图可视化
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

##### 信号网络的流形与分类学习 #####
# 回答问题：哪些信号通路在功能上或结构上是相似的
# 下载必要的Python umap-learn包
# 安装并加载 reticulate 包
# install.packages("reticulate")
# library(reticulate)
# py_install("umap-learn", pip = TRUE)

# 1. 功能相似性
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional") # 降维
cellchat <- netClustering(cellchat, type = "functional") # 聚类
netVisual_embedding(cellchat, type = "functional") # 画图

# 2. 结构相似性
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

#### 保存CellChat对象 ####
saveRDS(cellchat, file = "cellchat_single.rds")


