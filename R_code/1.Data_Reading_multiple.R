# 如果您觉得这个分享对您有所帮助，欢迎关注
# Bilibili：Broca区想发言
# GitHub：A-Ben-Mao

# GEO数据下载
# 官网：https://www.ncbi.nlm.nih.gov/geo/

# 本次所使用的数据集为
# GSE302285
# GSE291575
# GSE153935
# GSE165722

# 加载R包
# install.packages('Seurat')
library(Seurat)
library(data.table)
library(stringr)
library(tibble)

#### 1.matrix.mtx、genes.tsv和barcodes.tsv ####
# 最常见的数据格式
# 本示例使用GSE302285

# 需要读取的文件目录
setwd("文件夹目录")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files()
samples

# 创建一个空的列表
seurat_list <- list()

# 读取数据并创建Seurat对象
# 删除，小于200个基因表达的细胞，小于3个细胞表达的基因
for (sample in samples) {
  
  # 读取10x数据
  seurat_data <- Read10X(data.dir = sample)
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = seurat_data, # 原始表达矩阵
    project = sample,     # 为对象设置项目名，metadata中可区分来源
    min.features = 200,   # 细胞阈值
    min.cells = 3         # 基因阈值
  )
  
  # 添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

# 合并多个 Seurat 对象（跨样本合并）
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

# 将Layers融合
sc_data = JoinLayers(seurat_combined)

# 保存Seurat对象
saveRDS(sc_data, file = "Seurat_combined.rds")

#### 2.H5格式 #####
# 本示例使用GSE291575

# 需要读取的文件目录
setwd("文件夹目录")

# 列出当前目录下所有文件后缀为".h5"的文件
samples_h5=list.files(pattern = "\\.h5$")
samples_h5

# 依次读入并创建Seurat对象
sceList = lapply(samples_h5, function(x){
  a=Read10X_h5(x)
  p=str_split(x,'_',simplify = T)[,1]
  sce <- CreateSeuratObject(a,
                            project = p ,
                            min.features = 200,
                            min.cells = 3)
})

# 获取样本对应的GSM号
folders = substr(samples_h5,1,10)
folders

# 使用merge函数进行合并
sce.big <- merge(sceList[[1]], 
                 y = sceList[-1], 
                 add.cell.ids = folders)

# 将Layers融合
sc_data = JoinLayers(sce.big)

# 保存Seurat对象
saveRDS(sc_data, file = "Seurat_combined_h5.rds")

#### 以下内容仅经过整理后运行无报错，未进一步完善 ####
#### 3.TXT或CSV格式 ####
# 本示例使用GSE153935

# 需要读取的文件目录
setwd("文件夹目录")

# 使用fread函数读取数据，该函数读取速度较快
# 且这一类数据大多没有固定格式，每一个上传的文件都可能不一样
sc_data <-fread("GSE153935_TLDS_AllCells.txt", sep="\t")
sc_data[1:5,1:5] # 查看前五行五列数据

# 将第一列转换为行名
sc_data <-  column_to_rownames(sc_data,"V1")

# 创建seurat对象
sc_data <- CreateSeuratObject(sc_data, min.features = 300, min.cells = 3)

# 将Layers融合
sc_data = JoinLayers(sc_data)

# 保存Seurat对象
saveRDS(sc_data, file = "Seurat_combined_txt.rds")

#### 4.TSV格式 ####
# 本示例使用GSE165722

# 需要读取的文件目录
setwd("文件夹目录")

# 获取样本文件名
cellname=list.files(pattern = '.txt')
cellname
counts=list.files(pattern = '.tsv')
counts

# 获取样本名称
folders = substr(cellname,1,10)
folders

# 读取数据并创建Seurat对象
sceList = list()
for (i in 1:length(cellname)) {
  #读入
  abc123 = as.data.frame(fread(cellname[i]))
  abc456 = as.data.frame(fread(counts[i]))
  #将gene列转化为行名，基因名
  abc456 <-  column_to_rownames(abc456,"gene")
  #添加矩阵的列名，细胞名
  colnames(abc456) = abc123[,1]
  #创建seurat对象
  sce <- CreateSeuratObject(abc456,project = folders[i],min.features = 300, min.cells = 3)
  #依次放入list中
  sceList[i] = sce
}

# 使用merge函数进行合并
sce.big <- merge(sceList[[1]], 
                 y = sceList[-1], 
                 add.cell.ids = folders)

# 将Layers融合
sc_data = JoinLayers(sce.big)

# 保存Seurat对象
saveRDS(sc_data, file = "Seurat_combined_tsv.rds")

#### 5.RDS/RDATA文件 ####
# 需要读取的文件目录
setwd("文件夹目录")

# 读取RDA文件
load(file ="文件名")

# 读取RDS文件
sc_data = readRDS("文件名")

# 保存Seurat对象

saveRDS(sc_data, file = "Seurat_combined_r.rds")
