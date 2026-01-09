library(Biostrings)
library(msa)

setwd("/home/Dong/Global_COVID/genomic_analysis")
# 读取少量序列进行快速测试
fasta <- readDNAStringSet("gisaid_hcov-19_2025_07_18_12.fasta")

# 随机选择20条序列做快速测试
test_seq <- fasta[sample(length(fasta), 20)]

cat("开始速度测试...\n")
cat("测试序列数量:", length(test_seq), "\n")

# 定义safe_msa函数（简化版）
safe_msa <- function(sequences, max_seq = 500) {
  # MAFFT方法
  tryCatch({
    cat("测试MAFFT...\n")
    temp_file <- tempfile(fileext = ".fasta")
    writeXStringSet(sequences, temp_file)
    
    mafft_cmd <- paste("mafft --auto --quiet", temp_file)
    result_file <- tempfile(fileext = ".aln")
    
    start_time <- Sys.time()
    exit_code <- system(paste(mafft_cmd, ">", result_file))
    end_time <- Sys.time()
    
    if(exit_code == 0) {
      # 读取比对结果
      aligned_seq <- readDNAStringSet(result_file)
      unlink(c(temp_file, result_file))
      
      time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
      cat("MAFFT成功! 用时:", round(time_taken, 3), "秒\n")
      cat("MAFFT结果: ", length(aligned_seq), "条序列，长度", width(aligned_seq)[1], "\n")
      
      # 返回实际的比对结果和方法信息
      result <- list(
        alignment = aligned_seq,
        method = "MAFFT",
        time = time_taken
      )
      return(result)
    } else {
      stop("MAFFT失败")
    }
  }, error = function(e) {
    cat("MAFFT失败:", e$message, "\n")
    
    # 回退到ClustalOmega
    tryCatch({
      cat("测试ClustalOmega...\n")
      start_time <- Sys.time()
      res_align <- msa(sequences, method="ClustalOmega", type="dna")
      end_time <- Sys.time()
      time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
      cat("ClustalOmega成功! 用时:", round(time_taken, 3), "秒\n")
      cat("ClustalOmega结果: ", nrow(res_align), "条序列，长度", ncol(res_align), "\n")
      
      # 返回实际的比对结果和方法信息
      result <- list(
        alignment = res_align,
        method = "ClustalOmega", 
        time = time_taken
      )
      return(result)
    }, error = function(e2) {
      cat("ClustalOmega失败:", e2$message, "\n")
      return(list(
        alignment = NULL,
        method = "Failed",
        time = NA
      ))
    })
  })
}

# 运行测试
result <- safe_msa(test_seq)
cat("使用方法:", result$method, "\n")
cat("用时:", result$time, "秒\n")

# 显示比对结果的详细信息
if(!is.null(result$alignment)) {
  if(result$method == "MAFFT") {
    cat("MAFFT比对结果详情:\n")
    cat("- 序列数量:", length(result$alignment), "\n")
    cat("- 比对长度:", width(result$alignment)[1], "\n")
    cat("- 前3条序列名:", paste(names(result$alignment)[1:3], collapse=", "), "\n")
  } else if(result$method == "ClustalOmega") {
    cat("ClustalOmega比对结果详情:\n") 
    cat("- 序列数量:", nrow(result$alignment), "\n")
    cat("- 比对长度:", ncol(result$alignment), "\n")
  }
} else {
  cat("比对失败，无结果\n")
} 

# 访问实际的比对结果
cat("\n=== 比对结果对象 ===\n")
if(result$method == "MAFFT") {
  cat("MAFFT比对结果对象类型:", class(result$alignment), "\n")
  # 显示前几条序列的开头
  if(length(result$alignment) > 0) {
    cat("第一条序列前50个碱基:\n")
    cat(substr(as.character(result$alignment[[1]]), 1, 50), "\n")
  }
} else if(result$method == "ClustalOmega") {
  cat("ClustalOmega比对结果对象类型:", class(result$alignment), "\n")
}