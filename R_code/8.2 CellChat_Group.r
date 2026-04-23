# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面：
# 该数据大多是通过数据计算而得到对应的结果
# 但生物学意义远大于数量！

# 加载R包
library(CellChat)
library(patchwork)

# 设置工作目录
setwd("文件目录")

# 读取2个已独立完成分析的CellChat对象
cellchat.Control <- readRDS("对应文件地址")
cellchat.Disease <- readRDS("对应文件地址")

# 将多个对象放入列表，命名分组
object.list <- list(Control = cellchat.Control, Disease = cellchat.Disease)

# 合并多个CellChat对象
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 打印合并后的对象信息
cellchat

# 释放内存
rm(cellchat.Control,cellchat.Disease)
gc()

#### Part 1. 识别改变的相互作用和细胞群 ####
# 找细胞

##### 比较交互作用总数和交互作用强度 #####
# 宏观上快速判断两组样本的细胞通讯是整体增强还是减弱
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

##### 比较不同细胞群之间的相互作用数量和相互作用强度 ####
# 具体观察哪些细胞类型之间的通讯发生了显著变化

# (A)差异圈图：可视化两组细胞通讯的增减变化
# 结果为第二个数据集与第一个数据集的比较结果
# 红色为increased，蓝色为decreased
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# (B)差异热图：定量、精细化展示所有细胞对的差异
# 相比于圈图更加精准、定量
# 行 = 信号发送细胞(Outgoing)
# 列 = 信号接收细胞(Incoming)
# 色块颜色：红 = 增强，蓝 = 减弱；颜色深浅 = 变化幅度
# 顶部柱子 = 接收信号的总变化(哪类细胞接收的信号变化最大)
# 右侧柱子 = 发送信号的总变化(哪类细胞发送的信号变化最大)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# (C)独立圈图：每个分组单独画网络(统一尺度)
# 不对比差异，而是分别画出 Control 组、Disease 组的完整通讯圈图
# 用相同尺度直观对比
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# (D)聚合圈图：简化细胞类型，分析大类细胞通讯
# 如果细胞类型太多，网络太复杂，则可以考虑将细胞聚合为大类后再分析
# 手动定义细胞大类
group.cellType <- c(rep("A_cell", 3), rep("B_cell", 3), rep("C_cell", 3))
group.cellType <- factor(group.cellType, levels = c("A_cell", "B_cell", "C_cell"))
# 聚合细胞通讯
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
# 重新合并对象
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 绘制 Control / Disease 组的聚合通讯圈图
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# 绘制聚合后的差异圈图
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

##### 二维空间对比主要信号源和信号靶 #####
# 通过对比 传出(Outgoing)和传入(Incoming) 信号强度
# 快速识别：不同分组间，发送/接收信号发生显著变化的细胞群

# (A)识别发送或接收信号发生显著变化的细胞群
# 看图简易说明：
# 找位置变化：对照组在左下角(弱发+弱收)→疾病组跑到右上角(强发+强收)=关键驱动细胞
# 找大小变化：疾病组点明显变大 = 通讯全面增强

# 右上角：核心驱动细胞(既发信号又收信号)
# 左上角：纯信号源(只发不收)
# 右下角：纯信号靶(只收不发)
# 左下角：配角细胞(不参与核心通讯)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

# (B)识别特定细胞群的信号变化
# 在(A)中找到了关键变化的细胞，或有感兴趣的细胞
# 深入看特定细胞哪些信号通路的发送/接收发生了改变
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OPCs", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inhibitory_Neurons", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))

#### Part 2. 识别具有不同网络结构和相互作用强度的异常信号传导 ####
# 找信号通路

##### 识别显著差异的信号网络以及信号组 #####
# 识别网络结构/互作强度发生显著改变(Up Or Down)的信号通路
# 并根据功能相似性、结构相似性对信号通路分组
# 通过量化两组间信号网络的差异，找到疾病驱动的关键信号通路

# 加载环境
library(reticulate)
py_install("umap-learn", pip = TRUE)

# 基于功能相似性分析信号通路
# 看图简易说明：
# 每一个点 = 一个信号通路
# 同一个通路：在 Control、Disease 组各有一个位置
# 两点距离越远 = 该通路在两组间差异越大、重构越剧烈
# 点聚在一起 = 功能相似的通路
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional") # 计算功能相似性
cellchat <- netEmbedding(cellchat, type = "functional") # 降维
cellchat <- netClustering(cellchat, type = "functional") # 聚类
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# 基于结构相似性分析信号通路(可选)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural") # 降维
cellchat <- netClustering(cellchat, type = "structural") # 聚类
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

# 计算并可视化学习到的联合流形中的路径距离
# 每个信号通路在两组间的差异大小
# 距离越大，网络差异越大
# 仅计算两组共有的通路
rankSimilarity(cellchat, type = "functional")
rankSimilarity(cellchat, type = "structural")

##### 识别具有不同相互作用强度的改变的信号传导 #####
# 通过对比信息流(信号总强度)，识别两组间发生显著性改变的信号通路
# 关闭 (turn off)：对照组有信号，疾病组完全消失
# 减弱 (decrease)：疾病组信号强度显著下降
# 开启 (turn on)：对照组无信号，疾病组新出现
# 增强 (increase)：疾病组信号强度显著上升

# (A)对比信号通路的整体信息流
# 堆叠柱状图（两组信号强度叠加展示）
gg1 <- rankNet(
  cellchat, 
  mode = "comparison",  # 关键：分组对比模式
  measure = "weight",   # 关键：用信号强度（weight），不用数量
  stacked = T,          # 堆叠模式
  do.stat = TRUE        # 关键：做配对Wilcoxon检验，筛选显著差异通路
)
# 非堆叠柱状图（两组信号分开展示）
gg2 <- rankNet(
  cellchat, 
  mode = "comparison",
  measure = "weight",
  stacked = F,          # 非堆叠模式
  do.stat = TRUE
)
gg1 + gg2

# (B)对比细胞群的传出/传入信号模式
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
# 热图1：Outgoing 传出信号（细胞发送信号）
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = "Control", width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = "Disease", width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # 并排拼图
# 热图2：Incoming 传入信号（细胞接收信号）
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = "Control", width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = "Disease", width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# 热图3：All 全部信号（发送+接收）
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = "Control", width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = "Disease", width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#### Part 3. 识别上调和下调的信号配体-受体对 ####
# 找配体-受体对

##### 通过比较通讯概率识别异常信号 #####
# 查看细胞类型序号
cell_types <- unique(cellchat@idents[[2]])
cell_types <- sort(cell_types)
data.frame(序号 = 1:length(cell_types), 细胞类型 = cell_types)

# 全局对比气泡图
netVisual_bubble(
  cellchat, 
  sources.use = 3,       # 信号发送细胞：第X种细胞类型
  targets.use = c(1:9),  # 信号接收细胞：第a~b种细胞类型，可直接写名称，英文逗号隔开
  comparison = c(1, 2),  # 对比分组：1=Control, 2=Disease
  angle.x = 45           # X轴文字倾斜45度，防止重叠
)

# 分开展示上调 / 下调 LR 对
gg1 <- netVisual_bubble(cellchat, 
                        sources.use = 3, 
                        targets.use = c(1:9),  
                        comparison = c(1, 2), 
                        max.dataset = 2, 
                        title.name = "Increased signaling in Disease", 
                        angle.x = 45, 
                        remove.isolate = T
)
gg2 <- netVisual_bubble(cellchat, 
                        sources.use = 3, 
                        targets.use = c(1:9),  
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = "Decreased signaling in Disease", 
                        angle.x = 45, 
                        remove.isolate = T
)
gg1 + gg2

# 提取LR对
signaling_up = gg1$data
signaling_down = gg2$data

##### 利用差异表达分析识别功能异常的信号传导 #####
# 定义疾病组（Disease）
pos.dataset = "Disease"
# 定义结果存储名称
features.name = paste0(pos.dataset, ".merged")

# 鉴定差异表达基因（DEG）
cellchat <- identifyOverExpressedGenes(
  cellchat,
  group.dataset = "datasets",  # 分组依据
  pos.dataset = pos.dataset,   # 阳性组（疾病）
  features.name = features.name,
  thresh.pc = 0.1,      # 表达比例阈值
  thresh.fc = 0.05,     # 差异倍数阈值
  thresh.p = 0.05       # 显著性阈值
)

# 把DEG结果映射到已推断的通讯网络
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)

# 筛选上调/下调的 LR 对
net.up <- subsetCommunication(
  cellchat, net = net, datasets = "Disease",
  ligand.logFC = 0.05  # 配体logFC阈值
)
net.down <- subsetCommunication(
  cellchat, net = net, datasets = "Control",
  ligand.logFC = -0.05
)

# 提取上调信号的所有基因
# gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
# 提取下调信号的所有基因
# gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

##### 可视化已识别的上调和下调信号配体-受体对 #####
# 对通过DEG的分析方法得到的结果进行可视化

# (A)气泡图
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# (B)弦图
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 3, targets.use = c(1:9), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 3, targets.use = c(1:9), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# (C)词云
# install.packages("wordcloud")
library("wordcloud")
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)

#### Part 4. 可视化 ####
# 使用层级图、圆图或弦图直观地比较细胞间通讯

# 圈图对比
pathways.show <- c("PTN")  # 选择你要展示的信号通路
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]], 
    signaling = pathways.show,  # 指定通路
    layout = "circle",          # 圈图布局
    edge.weight.max = weight.max[1], # 统一最大边权重
    edge.width.max = 10, 
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

# 热图对比
pathways.show <- c("PTN")  # 选择你要展示的信号通路
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(
    object.list[[i]], 
    signaling = pathways.show, 
    color.heatmap = "Reds",
    title.name = paste(pathways.show, "signaling ",names(object.list)[i])
  )
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# 基础弦图对比
pathways.show <- c("PTN")  # 选择你要展示的信号通路
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]], 
    signaling = pathways.show, 
    layout = "chord",  # 弦图布局
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

# 细胞大类聚合弦图
group.cellType <- c(rep("A_cell", 3), rep("B_cell", 3), rep("C_cell", 3))
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("PTN")  # 选择你要展示的信号通路
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(
    object.list[[i]], 
    signaling = pathways.show, 
    group = group.cellType,  # 应用细胞大类
    title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i])
  )
}

# 指定发送细胞-接收细胞的精细可视化
# 下为几种细胞选择的写法
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(
    object.list[[i]], 
    sources.use = 3,        # 【发送】仅第3个细胞
    targets.use = c(1:9),   # 【接收】第1-9个细胞
    lab.cex = 0.5
  )}

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(
    object.list[[i]], 
    sources.use = c(1,2,3,4),  # 【发送】第1-4个细胞
    targets.use = c(8,10),     # 【接收】第8、10个细胞
    legend.pos.x = 10
  )}

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(
    object.list[[i]], 
    sources.use = c(1,2,3,4),  # 【发送】第1-4个细胞
    targets.use = c(5:9),      # 【接收】第5-9个细胞
    slot.name = "netP"         # 【修改】展示信号通路，不是LR对
  )}

#### Part 5. 比较不同数据集之间的信号基因表达分布 ####
# 绘制与 LR 对或信号通路相关的信号基因的基因表达分布
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Control", "Disease"))
plotGeneExpression(
  cellchat, 
  signaling = "PTN",    # 指定要画的信号通路
  split.by = "datasets",# 按分组拆分对比（Control/Disease）
  colors.ggplot = T,    # 使用ggplot2默认配色
  type = "violin"       # 绘图类型：小提琴图（展示表达分布）
)

#### 保存CellChat对象 ####
save(object.list, file = "cellchat_object.list__Ctrl_Disease.RData")
save(cellchat, file = "cellchat_merged_Ctrl_Disease.RData")

