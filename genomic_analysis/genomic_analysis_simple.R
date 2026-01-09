library(ape)
library(msa)
library(pegas)
library(devtools)
library(PopGenome)
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# 尝试加载可选包，如果失败则使用基础功能
phangorn_available <- requireNamespace("phangorn", quietly = TRUE)
if(phangorn_available) {
  library(phangorn)
  cat("✓ phangorn 包可用，将使用高级系统发育分析\n")
} else {
  cat("⚠️  phangorn 包不可用，将使用基础系统发育分析\n")
}

# 智能设置并行核心数
max_cores <- detectCores()
ncores <- min(8, max(1, max_cores - 2))
cat("检测到", max_cores, "个CPU核心，使用", ncores, "个核心进行并行计算\n")

# 尝试创建并行集群，如果失败则使用顺序处理
use_parallel <- TRUE
tryCatch({
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  cat("并行集群创建成功\n")
}, error = function(e) {
  cat("并行集群创建失败:", e$message, "\n将使用顺序处理\n")
  use_parallel <<- FALSE
  registerDoSEQ()
})

setwd("/home/Dong/Global_COVID/genomic_analysis")

fasta1 <- readDNAStringSet("gisaid_hcov-19_2025_07_18_12.fasta")
fasta2 <- readDNAStringSet("gisaid_hcov-19_2025_07_18_13.fasta")
fasta <- c(fasta1, fasta2)

# 添加序列长度检查
cat("总序列数量:", length(fasta), "\n")
cat("序列长度范围:", min(width(fasta)), "-", max(width(fasta)), "\n")

# 内存使用估算
approx_mem_gb <- (length(fasta) * mean(width(fasta)) * 8) / (1024^3)
cat("预估内存需求:", round(approx_mem_gb, 2), "GB\n")

# 提取时间信息
headers <- names(fasta)
dates <- sub(".*\\|([0-9]{4}-[0-9]{2}-[0-9]{2})$", "\\1", headers)
months <- substr(dates, 1, 7)
all_months <- sort(unique(months))
cat("发现", length(all_months), "个月份的数据\n")

# 简化的多序列比对函数 - 增强版
safe_msa <- function(sequences, max_seq = 100, method_preference = c("mafft", "clustalo", "muscle")) {
  if(length(sequences) > max_seq) {
    cat("序列数量过多 (", length(sequences), "), 随机采样", max_seq, "条序列\n")
    sequences <- sequences[sample(length(sequences), max_seq)]
  } 
  
  # 尝试使用MAFFT (如果系统有安装)
  if("mafft" %in% method_preference) {
    tryCatch({
      cat("尝试使用MAFFT进行多序列比对...\n")
      # 检查MAFFT是否可用
      mafft_available <- system("which mafft", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
      if(mafft_available) {
        # 创建临时文件
        temp_input <- tempfile(fileext = ".fasta")
        temp_output <- tempfile(fileext = ".fasta")
        
        # 写入序列
        writeXStringSet(sequences, temp_input)
        
                 # 运行MAFFT
         cmd <- paste("mafft --auto --quiet", temp_input, ">", temp_output)
        system(cmd)
        
        # 读取结果
        if(file.exists(temp_output) && file.size(temp_output) > 0) {
          aligned_seqs <- readDNAStringSet(temp_output)
          
          # 清理临时文件
          unlink(c(temp_input, temp_output))
          
          cat("✓ MAFFT比对成功\n")
          return(list(
            alignment = aligned_seqs,
            method = "MAFFT",
            type = "external"
          ))
        }
      }
    }, error = function(e) {
      cat("MAFFT失败:", e$message, "\n")
    })
  }
  
  # 降级到ClustalOmega
  tryCatch({
    cat("使用ClustalOmega进行多序列比对...\n")
    res_align <- msa(sequences, method="ClustalOmega")
    return(list(
      alignment = res_align,
      method = "ClustalOmega", 
      type = "msa_package"
    ))
  }, error = function(e) {
    cat("ClustalOmega失败，尝试MUSCLE:", e$message, "\n")
    tryCatch({
      res_align <- msa(sequences, method="Muscle")
      return(list(
        alignment = res_align,
        method = "MUSCLE",
        type = "msa_package"
      ))
    }, error = function(e2) {
      cat("MUSCLE也失败:", e2$message, "\n")
      return(NULL)
    })
  })
}

# 序列修剪函数 (模拟trimAL功能)
trim_alignment <- function(alignment, gap_threshold = 0.5) {
  cat("进行序列修剪 (gap阈值:", gap_threshold, ")...\n")
  
  tryCatch({
    if(is.list(alignment) && "alignment" %in% names(alignment)) {
      # 处理MAFFT结果
      if(alignment$type == "external") {
        align_matrix <- as.matrix(alignment$alignment)
      } else {
        # 处理msa包结果
        align_matrix <- as.matrix(msaConvert(alignment$alignment, type="seqinr::alignment"))
      }
    } else {
      # 直接的对齐结果
      align_matrix <- as.matrix(msaConvert(alignment, type="seqinr::alignment"))
    }
    
    # 计算每列的gap比例
    n_seqs <- nrow(align_matrix)
    gap_props <- apply(align_matrix, 2, function(col) {
      sum(col == "-" | col == "N" | is.na(col)) / n_seqs
    })
    
    # 保留gap比例低于阈值的列
    keep_cols <- gap_props <= gap_threshold
    n_removed <- sum(!keep_cols)
    
    if(sum(keep_cols) < 10) {
      cat("警告：修剪后序列过短，跳过修剪步骤\n")
      return(alignment)
    }
    
    trimmed_matrix <- align_matrix[, keep_cols, drop = FALSE]
    cat("序列修剪完成：移除", n_removed, "列，保留", sum(keep_cols), "列\n")
    
    # 转换回DNAStringSet
    trimmed_seqs <- DNAStringSet(apply(trimmed_matrix, 1, paste, collapse=""))
    names(trimmed_seqs) <- rownames(align_matrix)
    
    return(list(
      alignment = trimmed_seqs,
      method = paste(alignment$method, "+ trimAL-like"),
      type = "trimmed",
      original_length = ncol(align_matrix),
      trimmed_length = ncol(trimmed_matrix),
      removed_columns = n_removed
    ))
    
  }, error = function(e) {
    cat("序列修剪失败:", e$message, "，使用原始比对\n")
    return(alignment)
  })
}

# 提取时间信息函数
extract_sequence_dates <- function(seq_names) {
  cat("提取序列时间信息...\n")
  
  dates <- sapply(seq_names, function(name) {
    # 尝试多种日期格式
    date_patterns <- c(
      "[0-9]{4}-[0-9]{2}-[0-9]{2}",     # YYYY-MM-DD
      "[0-9]{4}/[0-9]{2}/[0-9]{2}",     # YYYY/MM/DD  
      "[0-9]{4}\\.[0-9]{2}\\.[0-9]{2}"  # YYYY.MM.DD
    )
    
    for(pattern in date_patterns) {
      match <- regexpr(pattern, name)
      if(match > 0) {
        date_str <- regmatches(name, match)
        # 标准化日期格式
        date_str <- gsub("[/\\.]", "-", date_str)
        return(date_str)
      }
    }
    return(NA)
  })
  
  valid_dates <- sum(!is.na(dates))
  cat("成功提取", valid_dates, "/", length(dates), "个序列的时间信息\n")
  
  return(dates)
}

# 增强的系统发育树构建函数
build_simple_phylo_tree <- function(sequences, sample_size = 100, method = "NJ", 
                                   use_trimming = TRUE, gap_threshold = 0.5) {
  cat("\n=== 构建系统发育树 ===\n")
  
  if(length(sequences) > sample_size) {
    cat("序列数量过多，采样", sample_size, "条序列\n")
    sequences <- sequences[sample(length(sequences), sample_size)]
  }
  
  # 提取时间信息
  seq_dates <- extract_sequence_dates(names(sequences))
  
  tryCatch({
    # 1. 多序列比对 (优先使用MAFFT)
    cat("进行多序列比对...\n")
    alignment_result <- safe_msa(sequences, max_seq = sample_size)
    if(is.null(alignment_result)) {
      cat("多序列比对失败\n")
      return(NULL)
    }
    
    # 2. 序列修剪 (可选)
    if(use_trimming) {
      alignment_result <- trim_alignment(alignment_result, gap_threshold)
    }
    
    # 3. 转换为DNAbin格式
    cat("转换为DNAbin格式...\n")
    if(alignment_result$type == "external" || alignment_result$type == "trimmed") {
      # 直接从DNAStringSet转换
      aligned_seqs <- alignment_result$alignment
      # 转换为字符矩阵再转DNAbin
      align_matrix <- as.matrix(aligned_seqs)
      dna_align <- as.DNAbin(align_matrix)
    } else {
      # 从msa对象转换
      dna_align <- as.DNAbin(msaConvert(alignment_result$alignment, type="seqinr::alignment"))
    }
    
    # 4. 计算距离矩阵
    cat("计算距离矩阵...\n")
    dist_matrix <- dist.dna(dna_align, model = "TN93", pairwise.deletion = TRUE)
    
    # 5. 构建系统发育树
    cat("构建系统发育树 (方法:", method, ")...\n")
    if(method == "NJ") {
      tree <- nj(dist_matrix)
    } else if(method == "UPGMA") {
      tree <- upgma(dist_matrix)
    } else if(method == "ML" && phangorn_available) {
      # 使用phangorn的最大似然法
      cat("使用最大似然法构建进化树...\n")
      phyDat_obj <- phyDat(dna_align)
      tree_init <- nj(dist_matrix)
      
      # 模型测试
      cat("进行模型选择...\n")
      models_to_test <- c("JC", "HKY", "GTR")
      best_aic <- Inf
      best_model <- "JC"
      
      for(model in models_to_test) {
        tryCatch({
          fit_test <- pml(tree_init, phyDat_obj, model = model)
          if(fit_test$logLik < best_aic) {
            best_aic <- fit_test$logLik  
            best_model <- model
          }
        }, error = function(e) {
          cat("模型", model, "测试失败\n")
        })
      }
      
      cat("选择模型:", best_model, "\n")
      fit <- pml(tree_init, phyDat_obj, model = best_model)
      fit <- optim.pml(fit, rearrangement = "stochastic", control = pml.control(trace = 0))
      tree <- fit$tree
      method <- paste("ML", best_model, sep = "-")
    } else {
      cat("使用NJ方法作为默认\n")
      tree <- nj(dist_matrix)
      method <- "NJ"
    }
    
    # 6. 中点定根
    if(!is.rooted(tree)) {
      tree <- midpoint(tree)
    }
    
    # 7. 计算分子钟信号 (简化版)
    molecular_clock_signal <- NA
    if(!all(is.na(seq_dates))) {
      tryCatch({
        valid_dates <- seq_dates[!is.na(seq_dates)]
        if(length(valid_dates) >= 10) {
          # 简单的分子钟检测：时间与距离的相关性
          date_nums <- as.numeric(as.Date(valid_dates))
          tip_names <- names(seq_dates)[!is.na(seq_dates)]
          
          if(length(tip_names) >= 10) {
            # 计算根到顶点距离
            root_to_tip_dist <- cophenetic(tree)[1, tip_names]
            
            if(length(root_to_tip_dist) == length(date_nums)) {
              correlation <- cor(date_nums, root_to_tip_dist, use = "complete.obs")
              molecular_clock_signal <- correlation
              cat("分子钟信号相关性:", round(correlation, 3), "\n")
            }
          }
        }
      }, error = function(e) {
        cat("分子钟信号计算失败:", e$message, "\n")
      })
    }
    
    cat("系统发育树构建完成! 包含", Ntip(tree), "个样本\n")
    
    return(list(
      tree = tree,
      alignment = dna_align,
      distance_matrix = dist_matrix,
      method = method,
      alignment_method = alignment_result$method,
      trimming_info = if(use_trimming && "removed_columns" %in% names(alignment_result)) {
        list(
          original_length = alignment_result$original_length,
          trimmed_length = alignment_result$trimmed_length,
          removed_columns = alignment_result$removed_columns
        )
      } else NULL,
      sequence_dates = seq_dates,
      molecular_clock_signal = molecular_clock_signal,
      analysis_info = list(
        n_sequences = length(sequences),
        n_valid_dates = sum(!is.na(seq_dates)),
        analysis_date = Sys.Date()
      )
    ))
    
  }, error = function(e) {
    cat("构建系统发育树失败:", e$message, "\n")
    return(NULL)
  })
}

# 增强的进化树绘制函数
plot_simple_phylo_tree <- function(phylo_result, prefix = "phylo", width = 12, height = 8, 
                                  include_time_annotation = TRUE) {
  if(is.null(phylo_result)) {
    cat("没有进化树数据可绘制\n")
    return(NULL)
  }
  
  cat("绘制系统发育树...\n")
  tree <- phylo_result$tree
  
  tryCatch({
    # 1. 基础系统发育树图
    pdf(paste0(prefix, "_tree.pdf"), width = width, height = height)
    plot(tree, main = paste("系统发育树 (", phylo_result$method, "方法)"),
         cex = 0.8, no.margin = TRUE)
    add.scale.bar()
    
    # 添加分析信息
    info_text <- paste(
      "对齐方法:", phylo_result$alignment_method,
      if(!is.null(phylo_result$trimming_info)) {
        paste("| 修剪:", phylo_result$trimming_info$removed_columns, "列移除")
      } else "",
      if(!is.na(phylo_result$molecular_clock_signal)) {
        paste("| 分子钟相关性:", round(phylo_result$molecular_clock_signal, 3))
      } else ""
    )
    mtext(info_text, side = 1, line = -2, cex = 0.7, adj = 0)
    dev.off()
    
    # 2. PNG格式 (高分辨率)
    png(paste0(prefix, "_tree.png"), width = width*100, height = height*100, res = 300)
    plot(tree, main = paste("系统发育树 (", phylo_result$method, "方法)"),
         cex = 0.8, no.margin = TRUE)
    add.scale.bar()
    mtext(info_text, side = 1, line = -2, cex = 0.7, adj = 0)
    dev.off()
    
    # 3. 圆形树图
    pdf(paste0(prefix, "_circular_tree.pdf"), width = 10, height = 10)
    plot(tree, type = "fan", main = "圆形系统发育树", cex = 0.6)
    add.scale.bar()
    dev.off()
    
    # 4. 时间标注树图 (如果有时间信息)
    if(include_time_annotation && !is.null(phylo_result$sequence_dates)) {
      dates <- phylo_result$sequence_dates
      valid_dates <- dates[!is.na(dates)]
      
      if(length(valid_dates) >= 10) {
        cat("绘制时间标注进化树...\n")
        
        # 创建颜色映射
        unique_months <- sort(unique(substr(valid_dates, 1, 7)))
        n_colors <- length(unique_months)
        colors <- rainbow(n_colors)
        names(colors) <- unique_months
        
        # 为每个序列分配颜色
        tip_colors <- rep("black", Ntip(tree))
        names(tip_colors) <- tree$tip.label
        
        for(i in 1:length(dates)) {
          if(!is.na(dates[i])) {
            month <- substr(dates[i], 1, 7)
            tip_colors[names(dates)[i]] <- colors[month]
          }
        }
        
        # 绘制时间标注树
        pdf(paste0(prefix, "_time_annotated_tree.pdf"), width = 14, height = 10)
        plot(tree, tip.color = tip_colors[tree$tip.label], 
             main = "时间标注系统发育树", cex = 0.8)
        add.scale.bar()
        
        # 添加图例
        legend("topright", legend = unique_months, fill = colors, 
               title = "采样月份", cex = 0.6, ncol = 2)
        dev.off()
        
        cat("时间标注树图已保存:", paste0(prefix, "_time_annotated_tree.pdf"), "\n")
      }
    }
    
    # 5. 保存Newick格式
    write.tree(tree, file = paste0(prefix, "_tree.newick"))
    
    # 6. 生成分析报告
    generate_phylo_report(phylo_result, prefix)
    
    cat("进化树图已保存:\n")
    cat("- ", prefix, "_tree.pdf/png\n")
    cat("- ", prefix, "_circular_tree.pdf\n") 
    cat("- ", prefix, "_tree.newick\n")
    cat("- ", prefix, "_analysis_report.txt\n")
    
    return(tree)
    
  }, error = function(e) {
    cat("绘制进化树失败:", e$message, "\n")
    return(NULL)
  })
}

# 生成系统发育分析报告
generate_phylo_report <- function(phylo_result, prefix = "phylo") {
  report_file <- paste0(prefix, "_analysis_report.txt")
  
  cat("生成分析报告:", report_file, "\n")
  
  report_content <- c(
    "=== 系统发育分析报告 ===",
    paste("分析日期:", Sys.Date()),
    paste("分析时间:", Sys.time()),
    "",
    "== 基本信息 ==",
    paste("序列数量:", phylo_result$analysis_info$n_sequences),
    paste("构建方法:", phylo_result$method),
    paste("对齐方法:", phylo_result$alignment_method),
    paste("进化树节点数:", Ntip(phylo_result$tree)),
    paste("内部节点数:", Nnode(phylo_result$tree)),
    paste("是否有根:", is.rooted(phylo_result$tree)),
    paste("树长度:", round(sum(phylo_result$tree$edge.length, na.rm = TRUE), 6)),
    "",
    "== 序列处理信息 ==",
    if(!is.null(phylo_result$trimming_info)) {
      c(
        paste("原始比对长度:", phylo_result$trimming_info$original_length, "bp"),
        paste("修剪后长度:", phylo_result$trimming_info$trimmed_length, "bp"),
        paste("移除列数:", phylo_result$trimming_info$removed_columns),
        paste("保留比例:", round(phylo_result$trimming_info$trimmed_length / phylo_result$trimming_info$original_length * 100, 1), "%")
      )
    } else {
      "未进行序列修剪"
    },
    "",
    "== 时间信息 ==",
    paste("包含时间信息的序列:", phylo_result$analysis_info$n_valid_dates),
    paste("时间信息完整性:", round(phylo_result$analysis_info$n_valid_dates / phylo_result$analysis_info$n_sequences * 100, 1), "%"),
    if(!is.na(phylo_result$molecular_clock_signal)) {
      paste("分子钟信号强度:", round(phylo_result$molecular_clock_signal, 4))
    } else {
      "分子钟信号: 无法计算"
    },
    "",
    "== 质量评估 ==",
    if(!is.na(phylo_result$molecular_clock_signal)) {
      paste("分子钟质量:",
            if(abs(phylo_result$molecular_clock_signal) > 0.5) "好 (>0.5)" 
            else if(abs(phylo_result$molecular_clock_signal) > 0.3) "中等 (0.3-0.5)"
            else "差 (<0.3)")
    } else "",
    paste("系统发育信息量:", 
          if(sum(phylo_result$tree$edge.length, na.rm = TRUE) > 0.01) "充足"
          else if(sum(phylo_result$tree$edge.length, na.rm = TRUE) > 0.001) "中等"
          else "较低"),
    "",
    "== 文件输出 ==",
    paste("-", prefix, "_tree.pdf: 系统发育树图 (PDF)"),
    paste("-", prefix, "_tree.png: 系统发育树图 (PNG)"),
    paste("-", prefix, "_circular_tree.pdf: 圆形系统发育树"),
    if(!is.null(phylo_result$sequence_dates) && phylo_result$analysis_info$n_valid_dates >= 10) {
      paste("-", prefix, "_time_annotated_tree.pdf: 时间标注系统发育树")
    } else "",
    paste("-", prefix, "_tree.newick: Newick格式树文件"),
    paste("-", prefix, "_analysis_report.txt: 本分析报告"),
    "",
    "== 建议 ==",
    if(!is.na(phylo_result$molecular_clock_signal) && abs(phylo_result$molecular_clock_signal) < 0.3) {
      "- 分子钟信号较弱，建议检查数据质量或增加序列数量"
    } else "",
    if(phylo_result$analysis_info$n_valid_dates / phylo_result$analysis_info$n_sequences < 0.8) {
      "- 时间信息不完整，建议完善序列的采样日期信息"
    } else "",
    if(Ntip(phylo_result$tree) < 20) {
      "- 样本数量较少，建议增加更多序列以提高分析可靠性"
    } else "",
    "",
    "=== 报告结束 ==="
  )
  
  # 过滤空行和NULL
  report_content <- report_content[!is.null(report_content) & report_content != ""]
  
  writeLines(report_content, report_file)
  
  cat("分析报告已保存:", report_file, "\n")
}

# 月度系统发育多样性分析
monthly_phylo_diversity <- function(fasta, months, all_months) {
  cat("\n=== 月度系统发育多样性分析 ===\n")
  
  phylo_diversity_stats <- data.frame()
  
  for(m in all_months) {
    month_idx <- which(months == m)
    if(length(month_idx) >= 10) {  # 至少10个样本
      
      # 限制序列数量
      if(length(month_idx) > 50) {
        month_idx <- sample(month_idx, 50)
      }
      
      month_seq <- fasta[month_idx]
      
      tryCatch({
        # 构建简单的NJ树计算系统发育多样性
        cat("处理月份", m, "- 序列数:", length(month_seq), "\n")
        
        month_align_result <- safe_msa(month_seq, max_seq = length(month_seq))
        if(!is.null(month_align_result)) {
          # 处理不同类型的比对结果
          if(month_align_result$type == "external" || month_align_result$type == "trimmed") {
            # 直接从DNAStringSet转换
            aligned_seqs <- month_align_result$alignment
            align_matrix <- as.matrix(aligned_seqs)
            month_dna <- as.DNAbin(align_matrix)
          } else {
            # 从msa对象转换
            month_dna <- as.DNAbin(msaConvert(month_align_result$alignment, type="seqinr::alignment"))
          }
          
          month_dist <- dist.dna(month_dna, model = "raw", pairwise.deletion = TRUE)
          month_tree <- nj(month_dist)
          
          # 计算系统发育多样性指标
          pd_value <- sum(month_tree$edge.length, na.rm = TRUE)
          mean_pairwise_dist <- mean(month_dist, na.rm = TRUE)
          
          phylo_stats <- data.frame(
            month = m,
            n_sequences = length(month_seq),
            phylogenetic_diversity = pd_value,
            mean_pairwise_distance = mean_pairwise_dist,
            tree_length = sum(month_tree$edge.length, na.rm = TRUE),
            alignment_method = month_align_result$method
          )
          
          phylo_diversity_stats <- rbind(phylo_diversity_stats, phylo_stats)
          cat("月份", m, "系统发育多样性:", round(pd_value, 4), "(", month_align_result$method, ")\n")
        }
      }, error = function(e) {
        cat("月份", m, "系统发育多样性计算失败:", e$message, "\n")
      })
    } else {
      cat("月份", m, "样本数不足 (", length(month_idx), " < 10)，跳过\n")
    }
  }
  
  return(phylo_diversity_stats)
}

################################################################################
# 执行分析
################################################################################

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("开始系统发育分析\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# 1. 全体样本分析（采样）
cat("\n### 全体样本系统发育分析 ###\n")
if(length(fasta) > 200) {
  # 分层采样
  sampled_indices <- c()
  month_counts <- table(months)
  
  for(m in names(month_counts)) {
    month_idx <- which(months == m)
    sample_size <- min(20, length(month_idx))
    if(sample_size >= 5) {
      sampled_idx <- sample(month_idx, sample_size)
      sampled_indices <- c(sampled_indices, sampled_idx)
    }
  }
  
  cat("从", length(fasta), "条序列中采样", length(sampled_indices), "条用于全体分析\n")
  global_sequences <- fasta[sampled_indices]
} else {
  global_sequences <- fasta
}

# 构建全体系统发育树 - 使用增强功能
global_phylo <- build_simple_phylo_tree(
  global_sequences, 
  sample_size = min(length(global_sequences), 100), 
  method = "ML",        # 使用最大似然法
  use_trimming = TRUE,  # 启用序列修剪
  gap_threshold = 0.5   # 50% gap阈值，符合参考标准
)

if(!is.null(global_phylo)) {
  # 绘制全体进化树
  plot_simple_phylo_tree(global_phylo, prefix = "global_phylo", include_time_annotation = TRUE)
  
  # 保存统计信息
  global_tree_stats <- data.frame(
    total_tips = Ntip(global_phylo$tree),
    tree_length = sum(global_phylo$tree$edge.length, na.rm = TRUE),
    method_used = global_phylo$method,
    alignment_method = global_phylo$alignment_method,
    molecular_clock_signal = global_phylo$molecular_clock_signal,
    is_rooted = is.rooted(global_phylo$tree),
    analysis_date = Sys.Date()
  )
  write.csv(global_tree_stats, "global_phylo_stats.csv", row.names = FALSE)
  
  cat("✅ 全体样本系统发育分析完成！\n")
  cat("分子钟信号强度:", 
      if(!is.na(global_phylo$molecular_clock_signal)) round(global_phylo$molecular_clock_signal, 3) else "无法计算", "\n")
}

# 2. 代表性月份分析 - 使用不同方法
cat("\n### 代表性月份系统发育分析 ###\n")
month_counts <- table(months)
representative_months <- names(month_counts[month_counts >= 30])
representative_months <- head(sort(representative_months), 5)

cat("选择", length(representative_months), "个代表性月份进行分析\n")

# 为不同月份使用不同方法进行对比
methods_to_test <- c("NJ", "ML", "UPGMA")
method_idx <- 1

for(month in representative_months) {
  cat("\n--- 分析月份:", month, "---\n")
  
  month_indices <- which(months == month)
  if(length(month_indices) > 30) {
    month_indices <- sample(month_indices, 30)
  }
  
  month_sequences <- fasta[month_indices]
  cat("该月份包含", length(month_sequences), "条序列\n")
  
  # 循环使用不同方法
  current_method <- methods_to_test[((method_idx - 1) %% length(methods_to_test)) + 1]
  cat("使用", current_method, "方法构建进化树\n")
  
  # 构建月份进化树
  monthly_phylo <- build_simple_phylo_tree(
    month_sequences, 
    sample_size = length(month_sequences),
    method = current_method,
    use_trimming = TRUE,
    gap_threshold = 0.5
  )
  
  if(!is.null(monthly_phylo)) {
    plot_simple_phylo_tree(
      monthly_phylo, 
      prefix = paste0("phylo_", month, "_", current_method),
      include_time_annotation = TRUE
    )
    cat("月份", month, "的系统发育分析完成 (", current_method, "方法)\n")
  }
  method_idx <- method_idx + 1
}

# 3. 系统发育多样性统计 - 使用增强版
phylo_diversity_stats <- monthly_phylo_diversity(fasta, months, all_months)

if(nrow(phylo_diversity_stats) > 0) {
  write.csv(phylo_diversity_stats, "monthly_phylogenetic_diversity.csv", row.names = FALSE)
  cat("系统发育多样性统计已保存\n")
  
  # 绘制趋势图 - 增强版
  tryCatch({
    # PDF版本
    pdf("phylogenetic_diversity_trend.pdf", width = 14, height = 8)
    par(mfrow = c(1, 2))
    
    # 系统发育多样性趋势
    plot(1:nrow(phylo_diversity_stats), phylo_diversity_stats$phylogenetic_diversity,
         type = "b", xlab = "月份序号", ylab = "系统发育多样性",
         main = "COVID-19系统发育多样性时间趋势", pch = 19, col = "blue")
    grid()
    
    # 平均成对距离趋势
    plot(1:nrow(phylo_diversity_stats), phylo_diversity_stats$mean_pairwise_distance,
         type = "b", xlab = "月份序号", ylab = "平均成对距离",
         main = "平均成对遗传距离时间趋势", pch = 19, col = "red")
    grid()
    
    dev.off()
    
    # PNG版本
    png("phylogenetic_diversity_trend.png", width = 1400, height = 800, res = 300)
    par(mfrow = c(1, 2))
    
    plot(1:nrow(phylo_diversity_stats), phylo_diversity_stats$phylogenetic_diversity,
         type = "b", xlab = "月份序号", ylab = "系统发育多样性", 
         main = "COVID-19系统发育多样性时间趋势", pch = 19, col = "blue")
    grid()
    
    plot(1:nrow(phylo_diversity_stats), phylo_diversity_stats$mean_pairwise_distance,
         type = "b", xlab = "月份序号", ylab = "平均成对距离",
         main = "平均成对遗传距离时间趋势", pch = 19, col = "red") 
    grid()
    
    dev.off()
    
    cat("系统发育多样性趋势图已保存\n")
  }, error = function(e) {
    cat("绘制趋势图失败:", e$message, "\n")
  })
  
  # 生成多样性报告
  diversity_summary <- data.frame(
    指标 = c("总月份数", "分析成功月份", "成功率", "平均系统发育多样性", "最大系统发育多样性", "最小系统发育多样性"),
    值 = c(
      length(all_months),
      nrow(phylo_diversity_stats),
      paste0(round(nrow(phylo_diversity_stats)/length(all_months)*100, 1), "%"),
      round(mean(phylo_diversity_stats$phylogenetic_diversity, na.rm = TRUE), 4),
      round(max(phylo_diversity_stats$phylogenetic_diversity, na.rm = TRUE), 4),
      round(min(phylo_diversity_stats$phylogenetic_diversity, na.rm = TRUE), 4)
    )
  )
  
  write.csv(diversity_summary, "phylogenetic_diversity_summary.csv", row.names = FALSE)
  cat("多样性分析摘要已保存\n")
}

# 安全关闭并行集群
if(use_parallel && exists("cl")) {
  tryCatch({
    stopCluster(cl)
    cat("并行集群已关闭\n")
  }, error = function(e) {
    cat("关闭并行集群时出错:", e$message, "\n")
  })
}

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("🎉 系统发育分析完成！\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\n📁 生成的文件包括:\n")
cat("📊 系统发育树文件:\n")
cat("  - global_phylo_tree.pdf/png: 全体系统发育树图\n")
cat("  - global_phylo_tree.newick: 全体系统发育树文件 (Newick格式)\n")
cat("  - global_phylo_time_annotated_tree.pdf: 时间标注进化树 (如果有时间信息)\n")
cat("  - phylo_*_tree.pdf/png: 月份系统发育树图\n")
cat("  - *_analysis_report.txt: 详细分析报告\n")

cat("\n📈 统计分析文件:\n")
cat("  - monthly_phylogenetic_diversity.csv: 月度系统发育多样性统计\n")
cat("  - phylogenetic_diversity_trend.pdf/png: 系统发育多样性趋势图\n")
cat("  - phylogenetic_diversity_summary.csv: 多样性分析摘要\n")
cat("  - global_phylo_stats.csv: 全体样本进化树统计\n")

cat("\n🔬 分析特色:\n")
cat("  ✅ 支持MAFFT高速比对 (如果系统已安装)\n")
cat("  ✅ 自动序列修剪 (trimAL-like, 50% gap阈值)\n")
cat("  ✅ 多种树构建方法 (NJ, ML, UPGMA)\n")
cat("  ✅ 自动模型选择 (ML方法)\n")
cat("  ✅ 分子钟信号检测\n")
cat("  ✅ 时间标注可视化\n")
cat("  ✅ 详细分析报告\n")

cat("\n💡 使用建议:\n")
cat("  1. 查看 *_analysis_report.txt 了解分析质量\n")
cat("  2. 使用 Newick 文件导入到 IQ-TREE, MEGA 等软件\n")
cat("  3. 如需更高级分析，请使用 IQ-TREE + TreeTime 流程\n")
cat("  4. 分子钟信号 >0.5 表示数据适合时间标定分析\n")

cat("\n🎯 与参考标准对比:\n")
cat("  MAFFT 比对:     ✅ 支持 (如果系统已安装)\n") 
cat("  trimAL 修剪:    ✅ 模拟实现 (50% gap 阈值)\n")
cat("  IQ-TREE ML:     ⚠️  基础ML (建议外部使用IQ-TREE)\n")
cat("  TempEst 检测:   ✅ 简化版分子钟检测\n")
cat("  TreeTime 标定:  ❌ 需外部工具\n") 