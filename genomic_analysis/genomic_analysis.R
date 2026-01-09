library(ape)
library(msa)
library(pegas)
library(devtools)
library(PopGenome)
library(Biostrings)
library(parallel)
library(doParallel)
library(foreach)

# 智能设置并行核心数
max_cores <- detectCores()
ncores <- max(8, max(1, max_cores - 2))
cat("检测到", max_cores, "个CPU核心，使用", ncores, "个核心进行并行计算\n")

# 尝试创建并行集群，如果失败则使用顺序处理
use_parallel <- TRUE
tryCatch({
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  cat("并行集群创建成功\n")
}, error = function(e) {
  cat("并行集群创建失败:", e$message, "\n将使用顺序处理\n")
  use_parallel <- FALSE
  registerDoSEQ()
})

setwd("/home/Dong/Global_COVID/genomic_analysis")

fasta1 <- readDNAStringSet("gisaid_hcov-19_2025_07_18_12.fasta")
fasta2 <- readDNAStringSet("gisaid_hcov-19_2025_07_18_13.fasta")
#fasta  <- readDNAStringSet("gisaid_hcov-19_2025_07_18_07.fasta")
fasta <- c(fasta1, fasta2)

#system('mafft gisaid_hcov-19_2025_07_18_07.fasta > aligned.fasta')


# 添加序列长度检查
cat("总序列数量:", length(fasta), "\n")
cat("序列长度范围:", min(width(fasta)), "-", max(width(fasta)), "\n")

# 内存使用估算
approx_mem_gb <- (length(fasta) * mean(width(fasta)) * 8) / (1024^3)
cat("预估内存需求:", round(approx_mem_gb, 2), "GB\n")

#res_align <- msa(fasta, method="ClustalW")  # 若要用mafft需正确配置,如下
#cmd <- paste("mafft --auto", "gisaid_hcov-19_2025_07_18_07.fasta", ">", "aligned.fasta")
#status <- system(cmd)
  

#res_align <- msa(fasta, method="Muscle")
#alignment <- as.DNAbin(as.matrix(res_align))
#seg_sites <- length(seg.sites(alignment))
#cat("多态位点数（segregating sites）:", seg_sites, "\n")
#haps <- haplotype(alignment)
#hd <- hap.div(haps)
#cat("单倍型多样性(haplotype diversity):", hd, "\n")
#pi <- nuc.div(alignment)
#cat("核苷酸多样性(nucleotide diversity, pi):", pi, "\n")
#save(list = ls(), file = "/home/Dong/Global_COVID/genomic_analysis/genomic_analysis.Rdata")
########################################################################################
########################################################################################
headers <- names(fasta)
dates <- sub(".*\\|([0-9]{4}-[0-9]{2}-[0-9]{2})$", "\\1", headers)
months <- substr(dates, 1, 7)
all_months <- sort(unique(months))
cat("发现", length(all_months), "个月份的数据\n")

# 添加安全的多序列比对函数
safe_msa <- function(sequences) {
  # 限制序列数量避免崩溃
  # 写入临时fasta文件
  tmp_in <- tempfile(fileext = ".fasta")
  tmp_out <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(sequences, tmp_in)
  
  # 调用mafft
  cmd <- paste("mafft --auto --quiet", tmp_in, ">", tmp_out)
  status <- system(cmd)
  
  if(status != 0) {
    cat("mafft运行失败，状态码：", status, "\n")
    unlink(c(tmp_in, tmp_out))
    return(NULL)
  }
    aligned <- Biostrings::readDNAStringSet(tmp_out)  
    # 清理临时文件
    unlink(c(tmp_in, tmp_out))
    cat("MAFFT比对完成，得到", length(aligned), "条比对序列\n")
    return(aligned)
}

# 处理每月多样性分析
cat("开始月度多样性分析...\n")
if(use_parallel) {
  # 并行处理
  div_list <- foreach(m = all_months, .combine=rbind, 
                      .packages=c("msa", "ape", "pegas", "Biostrings"),
                      .errorhandling = "remove",
                      .export = c("fasta", "months", "safe_msa")) %dopar% {    
    tryCatch({
      idx <- which(months == m)
      if(length(idx) < 2) return(NULL)      
      cat("处理月份:", m, ", 序列数量:", length(idx), "\n")     
      fasta_sub <- fasta[idx]
      res_align <- safe_msa(fasta_sub)     
      if(is.null(res_align)) {
        cat("月份", m, "比对失败，跳过\n")
        return(NULL)
      }     
      #alignment <- as.DNAbin(as.matrix(res_align))
      alignment <- as.DNAbin(as.matrix(res_align))
  seg_sites <- length(seg.sites(alignment))
  haps <- haplotype(alignment)
  n_haps <- length(haps)
  hd <- hap.div(haps)
  pi <- nuc.div(alignment)
  theta_w <- theta.s(alignment)
  
  # 距离相关
  dist_matrix <- dist.dna(alignment, model="raw", pairwise.deletion=TRUE)
  mean_dist <- mean(dist_matrix, na.rm=TRUE)
  max_dist <- max(dist_matrix, na.rm=TRUE)
  
  # Tajima's D
  tajima_result <- tajima.test(alignment)
  tajima_d <- tajima_result$D
  
  # 有效群体大小估计 (基于 π 和 θw)
  ne_pi <- pi / (4 * 2.3e-6)  # 假设突变率 2.3e-6 per site per year
  ne_theta <- theta_w / (4 * 2.3e-6)
      
  div_result <- data.frame(
    month = m,
    n_seq = length(idx),
    seg_sites = seg_sites,
    n_haplotypes = n_haps,
    hap_diversity = hd,
    nuc_diversity_pi = pi,
    theta_watterson = theta_w,
    mean_pairwise_dist = mean_dist,
    max_pairwise_dist = max_dist,
    tajima_d = tajima_d,
    ne_estimate_pi = ne_pi,
    ne_estimate_theta = ne_theta,
        stringsAsFactors=FALSE
      )
    }, error = function(e) {
      cat("月份", m, "分析出错:", e$message, "\n")
      return(NULL)
    })
  }
} else {
  # 顺序处理
  div_list <- data.frame()
  for(m in all_months) {
    tryCatch({
      idx <- which(months == m)
      if(length(idx) < 2) next
      
      cat("处理月份:", m, ", 序列数量:", length(idx), "\n")
      
      fasta_sub <- fasta[idx]
      res_align <- safe_msa(fasta_sub)
      
      if(is.null(res_align)) {
        cat("月份", m, "比对失败，跳过\n")
        next
      }
      
      alignment <- as.DNAbin(as.matrix(res_align))
  seg_sites <- length(seg.sites(alignment))
  haps <- haplotype(alignment)
  n_haps <- length(haps)
  hd <- hap.div(haps)
  pi <- nuc.div(alignment)
  theta_w <- theta.s(alignment)
  
  # 距离相关
  dist_matrix <- dist.dna(alignment, model="raw", pairwise.deletion=TRUE)
  mean_dist <- mean(dist_matrix, na.rm=TRUE)
  max_dist <- max(dist_matrix, na.rm=TRUE)
  
  # Tajima's D
  tajima_result <- tajima.test(alignment)
  tajima_d <- tajima_result$D
  
  # 有效群体大小估计 (基于 π 和 θw)
  ne_pi <- pi / (4 * 2.3e-6)  # 假设突变率 2.3e-6 per site per year
  ne_theta <- theta_w / (4 * 2.3e-6)
      
  div_result <- data.frame(
        month = m,
        n_seq = length(idx),
        seg_sites = seg_sites,
    n_haplotypes = n_haps,
    hap_diversity = hd,
    nuc_diversity_pi = pi,
    theta_watterson = theta_w,
    mean_pairwise_dist = mean_dist,
    max_pairwise_dist = max_dist,
    tajima_d = tajima_d,
    ne_estimate_pi = ne_pi,
    ne_estimate_theta = ne_theta,
        stringsAsFactors=FALSE
      )
      div_list <- rbind(div_list, div_result)
    }, error = function(e) {
      cat("月份", m, "分析出错:", e$message, "\n")
    })
  }
}
write.csv(div_list, "monthly_diversity_stats.csv", row.names=FALSE)
# 处理月度多样性结果
if(is.null(div_list)) {
  cat("警告：月度多样性分析没有返回任何结果\n")
  div_list <- data.frame()
} else {
  if(is.data.frame(div_list)) {
    cat("成功分析", nrow(div_list), "个月份的多样性数据\n")
  } else {
    cat("警告：多样性分析结果格式异常\n")
    div_list <- data.frame()
  }
}

cat("div_list:\n")
print(div_list)
#write.csv(div_list, "monthly_diversity_stats.csv", row.names=FALSE)

# 抗原漂移分析 (Antigenic Drift Analysis)
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("开始抗原漂移分析\n")
cat(paste(rep("=", 60), collapse=""), "\n")

if(use_parallel) {
  # 并行处理
  adrift_list <- foreach(m = all_months, .combine = rbind, 
                         .packages = c("msa", "ape", "Biostrings"),
                         .errorhandling = "remove",
                         .export = c("fasta", "months", "safe_msa", "all_months")) %dopar% {
    
    tryCatch({
      idx_i <- which(months == m)
      if(length(idx_i) < 1) return(NULL)
      prev_months <- all_months[which(all_months < m)]
      prev_months <- tail(prev_months, 3)
      idx_j <- which(months %in% prev_months)
      if(length(idx_j) < 1) return(NULL)
      
      
      fasta_i <- fasta[idx_i]
      fasta_j <- fasta[idx_j]
      
      combined <- c(fasta_i, fasta_j)
      res_align <- safe_msa(combined)
      
      if(is.null(res_align)) {
        cat("月份", m, "的抗原漂移分析失败\n")
        return(NULL)
      }
      
      alignment <- as.DNAbin(as.matrix(res_align))
      dist_matrix <- dist.dna(alignment, model="N", pairwise.deletion=TRUE)

      n_i <- length(idx_i)
      n_j <- length(idx_j)
      
      # 安全的矩阵访问
      if(n_i > 0 && n_j > 0) {
        dist_sub <- as.matrix(dist_matrix)[1:n_i, (n_i+1):(n_i+n_j)]
        mean_dist <- mean(dist_sub, na.rm = TRUE)
        median_dist <- median(dist_sub, na.rm = TRUE)
        
        data.frame(month = m, mean_distance = mean_dist, median_distance = median_dist, stringsAsFactors=FALSE)
      } 
    }, error = function(e) {
      cat("月份", m, "抗原漂移分析出错:", e$message, "\n")
      return(NULL)
    })
  }
} else {
  # 顺序处理
  adrift_list <- data.frame()
  for(m in all_months) {
    tryCatch({
      idx_i <- which(months == m)
      if(length(idx_i) < 1) next
      prev_months <- all_months[which(all_months < m)]
      prev_months <- tail(prev_months, 3)
      idx_j <- which(months %in% prev_months)
      if(length(idx_j) < 1) next
      
      # 限制序列数量
      if(length(idx_i) > 350) idx_i <- sample(idx_i, 350)
      if(length(idx_j) > 1050) idx_j <- sample(idx_j, 1050)
      
      fasta_i <- fasta[idx_i]
      fasta_j <- fasta[idx_j]
      
      combined <- c(fasta_i, fasta_j)
      res_align <- safe_msa(combined)
      
      if(is.null(res_align)) {
        cat("月份", m, "的抗原漂移分析失败\n")
        next
      }
      
      alignment <- as.DNAbin(as.matrix(res_align))
      dist_matrix <- dist.dna(alignment, model="N", pairwise.deletion=TRUE)
      n_i <- length(idx_i)
      n_j <- length(idx_j)
      
      # 安全的矩阵访问
      if(n_i > 0 && n_j > 0) {
        dist_sub <- as.matrix(dist_matrix)[1:n_i, (n_i+1):(n_i+n_j)]
        mean_dist <- mean(dist_sub, na.rm = TRUE)
        median_dist <- median(dist_sub, na.rm = TRUE)
        
        adrift_result <- data.frame(month = m, mean_distance = mean_dist, median_distance = median_dist, stringsAsFactors=FALSE)
        adrift_list <- rbind(adrift_list, adrift_result)
      }
    }, error = function(e) {
      cat("月份", m, "抗原漂移分析出错:", e$message, "\n")
    })
  }
}
write.csv(adrift_list, "monthly_antigenic_drift.csv", row.names=FALSE)
# 处理抗原漂移结果
if(is.null(adrift_list)) {
  cat("警告：抗原漂移分析没有返回任何结果\n")
  adrift_list <- data.frame()
} else {
  if(is.data.frame(adrift_list)) {
    cat("成功分析", nrow(adrift_list), "个月份的抗原漂移\n")
  } else {
    cat("警告：抗原漂移结果格式异常\n")
    adrift_list <- data.frame()
  }
}

cat("adrift_list:\n")
print(adrift_list)
#write.csv(adrift_list, "monthly_antigenic_drift.csv", row.names=FALSE)

# 检测显著抗原漂移事件
if(!is.null(adrift_list) && is.data.frame(adrift_list) && nrow(adrift_list) > 0) {
  drift_mean <- adrift_list$mean_distance
  drift_z <- scale(drift_mean)
  threshold <- 0.5
  
  important_events <- adrift_list[abs(drift_z) > threshold, ]
  if(nrow(important_events) > 0) {
    cat("检测到显著抗原漂移事件:\n")
    print(important_events)
    write.csv(important_events, "significant_antigenic_drift_events.csv", row.names=FALSE)
  } else {
    cat("未检测到显著抗原漂移事件\n")
  }
} else {
  cat("没有足够的数据进行抗原漂移事件检测\n")
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
write.csv(important_events, "significant_antigenic_drift_events.csv", row.names=FALSE)
cat("\n分析完成！\n")
cat("生成的文件:\n")
cat("- monthly_diversity_stats.csv: 月度多样性统计\n")
cat("- monthly_antigenic_drift.csv: 月度抗原漂移分析\n")
cat("- significant_antigenic_drift_events.csv: 显著抗原漂移事件\n")
















