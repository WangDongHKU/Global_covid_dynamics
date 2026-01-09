# 系统发育分析R包安装脚本
# Phylogenetic Analysis R Package Installation Script

cat("=== 安装系统发育分析所需的R包 ===\n")

# 安装CRAN包
cran_packages <- c(
  "phangorn"    # 系统发育重建
)

cat("安装CRAN包...\n")
for(pkg in cran_packages) {
  if(!require(pkg, character.only = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "已安装\n")
  }
}

# 安装Bioconductor包
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "ggtree",      # 进化树可视化
  "treeio"       # 进化树文件处理
)

cat("安装Bioconductor包...\n")
for(pkg in bioc_packages) {
  if(!require(pkg, character.only = TRUE)) {
    cat("安装", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE)
  } else {
    cat(pkg, "已安装\n")
  }
}

# 验证所有包是否成功安装
cat("\n=== 验证包安装 ===\n")
all_packages <- c(cran_packages, bioc_packages)
missing_packages <- c()

for(pkg in all_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat("✗", pkg, "安装失败\n")
  } else {
    cat("✓", pkg, "安装成功\n")
  }
}

if(length(missing_packages) == 0) {
  cat("\n🎉 所有包安装成功！\n")
  cat("现在可以运行系统发育分析了。\n")
} else {
  cat("\n⚠️  以下包安装失败:\n")
  for(pkg in missing_packages) {
    cat("-", pkg, "\n")
  }
  cat("\n请手动安装失败的包或检查网络连接。\n")
}

# 显示关键包的版本信息
cat("\n=== 关键包版本信息 ===\n")
key_packages <- c("ape", "phangorn", "ggtree", "Biostrings", "msa")
for(pkg in key_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(pkg, ":", as.character(version), "\n")
  }
}

cat("\n安装完成！\n") 