# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# 写在前面：为什么要做 GO 和 KEGG 富集分析？
# 在先前的步骤中，我们得到了一堆差异基因（Marker Genes）
# 而富集分析是为了把这些基因注释到具体的生物学通路上
# 1. GO (Gene Ontology)：
#    - BP (生物学过程): 细胞在干什么？(如：神经元突触传递、免疫反应)
#    - CC (细胞组分): 蛋白在哪里？(如：细胞核、细胞膜)
#    - MF (分子功能): 分子有什么功能？(如：ATP结合)
# 2. KEGG (Pathway):
#    - 具体的信号通路图 (如：Alzheimer disease pathway, Calcium signaling pathway)
#    - 帮助我们要研究具体的分子机制和疾病关联。

#### 1. 环境准备与包加载 ####
# 加载R包
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 设置工作目录
setwd("文件目录")

#### 2. 数据读取与细胞类型选择 ####
# 读入文件
markers <- read_csv("04.Excitatory_deg.csv")

# 规范化 cluster 名字（防止文件名出错）
markers$cluster <- str_replace_all(markers$cluster, "/", ".")

# 查看当前有哪些 Cluster (方便你复制名字)
print(unique(markers$cluster))

# 输入你想分析的 Cluster
target_cluster <- "Excitatory"   # 根据研究进行修改

#### 3. 提取基因与 ID 转换 ####
# 提取基因
input_gene_symbol <- markers %>% 
  filter(cluster == target_cluster) %>% 
  pull(gene)

# 转换 ID (Symbol -> Entrez ID)
gene_map <- bitr(input_gene_symbol, 
                 fromType = "SYMBOL", 
                 toType   = "ENTREZID", 
                 OrgDb    = org.Hs.eg.db)
gene_entrez <- gene_map$ENTREZID

#### 4. GO 富集分析 (Gene Ontology) ####
message("--- 正在进行 GO 分析 ---")

# 核心分析代码
ego <- enrichGO(gene          = gene_entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",      # "BP", "CC", "MF" 或 "ALL"
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)       # 结果中自动将 ID 转回 Symbol

# 保存原始结果
ego_res <- as.data.frame(ego)
write.table(ego_res, 
            file = paste0(target_cluster, "_GO.csv"),
            sep = ",",
            quote = FALSE, 
            row.names = FALSE)

# 筛选后并保存
# 可根据"pvalue"(原始P值)、"p.adjust"（矫正P值）和"qvalue"（Q值，类似FDR）进行筛选
ego_filtered <- ego_res %>% filter(pvalue < 0.05 & qvalue < 0.2)
write.table(ego_filtered, 
            file = paste0(target_cluster, "_GO_filtered.csv"),  # 修改文件扩展名
            sep = ",",
            quote = FALSE, 
            row.names = FALSE)

#### 5. GO 的可视化 ####
# 1. 柱状图
# pdf(file = paste0(target_cluster, ".GObarplot.pdf"), width = 12, height = 8)
print(barplot(ego, drop = TRUE, showCategory = 10, split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scale = 'free'))
# dev.off()

# 2. 气泡图
# pdf(file = paste0(target_cluster, ".GObubble.pdf"), width = 12, height = 8)
print(dotplot(ego, showCategory = 10, orderBy = "GeneRatio", split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scale = 'free'))
# dev.off()

#### 6. KEGG 富集分析 (Pathway) ####
message("--- 正在进行 KEGG 分析 ---")

# 需要联网分析，可能由于网络原因报错，可多尝试几次
# 下列代码可能在一定程度上解决 KEGG 网络问题
# R.utils::setOption("clusterProfiler.download.method",'auto')
options(timeout = 99999) # 延长网络连接时间

# 核心代码
ekegg <- enrichKEGG(gene         = gene_entrez,
                    organism     = "hsa", # 选择对应物种，小鼠需要切换为mmu
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

# ID 翻译回 Symbol
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 保存原始结果
kegg_res <- as.data.frame(ekegg)
write.table(kegg_res, 
            file = paste0(target_cluster, "_KEGG.csv"),
            sep = ",",
            quote = FALSE, 
            row.names = FALSE)

# 筛选后并保存
# 可根据"pvalue"(原始P值)、"p.adjust"（矫正P值）和"qvalue"（Q值，类似FDR）进行筛选
kegg_filtered <- kegg_res %>% filter(pvalue < 0.05 & qvalue < 0.2)
write.table(kegg_filtered, 
            file = paste0(target_cluster, "_KEGG_filtered.csv"),
            sep = ",",
            quote = FALSE, 
            row.names = FALSE)

#### 7. KEGG 的可视化 ####
# 1. 柱状图
# pdf(file = paste0(target_cluster, ".KEGGbarplot.pdf"), width = 10, height = 7)
print(barplot(ekegg, drop = TRUE, showCategory = 15))
# dev.off()

# 2. 气泡图
# pdf(file = paste0(target_cluster, ".KEGGbubble.pdf"), width = 10, height = 7)
print(dotplot(ekegg, showCategory = 15, orderBy = "GeneRatio"))
# dev.off()

