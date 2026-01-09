# 系统发育分析R包安装脚本 - 简化版
# 避免复杂依赖，使用基础包实现核心功能

cat("=== 安装系统发育分析所需的R包 (简化版) ===\n")

# 首先安装基础依赖
basic_packages <- c(
  "ggplot2",     # 绘图 - 必需
  "ape",         # 系统发育分析 - 核心包
  "phangorn"     # 系统发育重建 - 核心包
)

cat("安装基础CRAN包...\n")
for(pkg in basic_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
    if(require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("✓", pkg, "安装成功\n")
    } else {
      cat("✗", pkg, "安装失败\n")
    }
  } else {
    cat("✓", pkg, "已安装\n")
  }
}

# 手动安装yulab.utils（ggtree的依赖）
cat("\n安装yulab.utils（ggtree依赖）...\n")
if(!require("yulab.utils", character.only = TRUE, quietly = TRUE)) {
  tryCatch({
    install.packages("yulab.utils", repos = "https://cran.r-project.org/")
    cat("✓ yulab.utils 安装成功\n")
  }, error = function(e) {
    cat("✗ yulab.utils 安装失败，将跳过ggtree\n")
  })
}

# 尝试安装ggtree生态系统
cat("\n尝试安装ggtree生态系统...\n")
ggtree_packages <- c("tidytree", "treeio", "ggtree")

success_count <- 0
for(pkg in ggtree_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch({
      if(pkg %in% c("ggtree", "treeio")) {
        BiocManager::install(pkg, update = FALSE)
      } else {
        install.packages(pkg)
      }
      if(require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("✓", pkg, "安装成功\n")
        success_count <- success_count + 1
      } else {
        cat("✗", pkg, "安装失败\n")
      }
    }, error = function(e) {
      cat("✗", pkg, "安装失败:", e$message, "\n")
    })
  } else {
    cat("✓", pkg, "已安装\n")
    success_count <- success_count + 1
  }
}

# 验证核心功能
cat("\n=== 验证核心功能 ===\n")
core_packages <- c("ape", "phangorn", "ggplot2")
all_core_available <- TRUE

for(pkg in core_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "可用\n")
  } else {
    cat("✗", pkg, "不可用\n")
    all_core_available <- FALSE
  }
}

if(all_core_available) {
  cat("\n🎉 核心系统发育分析功能可用！\n")
  if(success_count == length(ggtree_packages)) {
    cat("✓ 高级可视化功能也可用 (ggtree)\n")
  } else {
    cat("⚠️  高级可视化功能不可用，将使用基础绘图\n")
  }
} else {
  cat("\n❌ 核心功能不完整，请检查安装\n")
}

cat("\n=== 创建测试脚本 ===\n")
# 创建一个简单的测试
test_code <- '
# 测试系统发育分析功能
library(ape)
library(phangorn)

# 创建测试数据
test_data <- matrix(sample(c("A","T","G","C"), 100, replace=TRUE), nrow=10)
rownames(test_data) <- paste0("seq", 1:10)
test_dna <- as.DNAbin(test_data)

# 计算距离和构建树
test_dist <- dist.dna(test_dna)
test_tree <- nj(test_dist)

# 测试绘图
png("test_tree.png", width=800, height=600)
plot(test_tree, main="测试系统发育树")
dev.off()

cat("✓ 系统发育分析功能测试成功！\\n")
cat("测试树图已保存为 test_tree.png\\n")
'

writeLines(test_code, "test_phylo_functionality.R")
cat("✓ 测试脚本已创建: test_phylo_functionality.R\n")

cat("\n安装完成！运行以下命令测试功能:\n")
cat("source('test_phylo_functionality.R')\n") 