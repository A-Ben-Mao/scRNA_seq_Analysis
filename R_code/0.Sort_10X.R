# 原始目录
raw_dir <- "文件夹目录"

# 获取所有文件
files <- list.files(raw_dir, full.names = FALSE)

# 遍历文件
for (file in files) {
  
  # 只处理 gz 文件
  if (!grepl("\\.gz$", file)) next
  
  # 完整路径
  file_path <- file.path(raw_dir, file)
  
  # 判断文件类型
  if (grepl("barcodes.tsv.gz$", file)) {
    file_type <- "barcodes.tsv.gz"
  } else if (grepl("features.tsv.gz$", file)) {
    file_type <- "features.tsv.gz"
  } else if (grepl("matrix.mtx.gz$", file)) {
    file_type <- "matrix.mtx.gz"
  } else {
    next
  }
  
  # 提取样本名（去掉后三种后缀）
  sample_name <- sub("_(barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz)$", "", file)
  
  # 样本文件夹（在 RAW 上一层目录）
  parent_dir <- dirname(raw_dir)
  sample_folder <- file.path(parent_dir, sample_name)
  
  # 如果不存在则创建
  if (!dir.exists(sample_folder)) {
    dir.create(sample_folder)
  }
  
  # 目标文件路径（复制并重命名）
  dest_path <- file.path(sample_folder, file_type)
  
  # 复制文件（保留原文件）
  file.copy(file_path, dest_path, overwrite = TRUE)
  
  cat("Copied:", file, "→", dest_path, "\n")
}

cat("文件整理完成！（原文件已保留）\n")

