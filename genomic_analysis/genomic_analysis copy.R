library(ape)
library(msa)
library(pegas)
library(devtools)
library(PopGenome)
library(Biostrings)
#BiocManager::install("vegan")
#BiocManager::install("adegenet")

#install.packages("poppr")
#library(poppr)
#install.packages("adegenet")
#library(adegenet)
#devtools::install_github("pievos101/PopGenome")


setwd("/home/Dong/Global_COVID/genomic_analysis")

fasta1 <- readDNAStringSet("gisaid_hcov-19_2025_07_18_12.fasta")
fasta2 <- readDNAStringSet("gisaid_hcov-19_2025_07_18_13.fasta")
fasta <- c(fasta1, fasta2)

#res_align <- msa(fasta, method="ClustalW")  # 若要用mafft需正确配置,如下
#res_align <- msa(fasta, method="Muscle")
#alignment <- as.DNAbin(msaConvert(res_align, type="seqinr::alignment"))
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

results <- data.frame(
  month = character(),
  n_seq = integer(),
  seg_sites = integer(),
  hap_div = numeric(),
  nuc_div = numeric(),
  stringsAsFactors = FALSE
)

all_months <- sort(unique(months))
for(m in all_months) {
  idx <- which(months == m)
  if(length(idx) < 2) next   # 少于2条序列跳过
  fasta_sub <- fasta[idx]
  
  # 多序列比对
  res_align <- msa(fasta_sub, method="Muscle")
  alignment <- as.DNAbin(msaConvert(res_align, type="seqinr::alignment"))
  seg_sites <- length(seg.sites(alignment))
  haps <- haplotype(alignment)
  hd <- hap.div(haps)
  pi <- nuc.div(alignment)
  results <- rbind(results,
    data.frame(
      month = m,
      n_seq = length(idx),
      seg_sites = seg_sites,
      hap_div = hd,
      nuc_div = pi))
  cat(m, ": n=", length(idx), " S=", seg_sites, " Hd=", hd, " pi=", pi, "\n")
}
write.csv(results, "monthly_diversity_stats.csv", row.names=FALSE)

# --- Antigenic drift analysis: pairwise distances between current month and previous 6 months ---
antigenic_drift_results <- data.frame(
  month = character(),
  mean_distance = numeric(),
  median_distance = numeric(),
  stringsAsFactors = FALSE
)

for(m in all_months) {
  idx_i <- which(months == m)
  if(length(idx_i) < 1) next
  # Get previous 6 months
  prev_months <- all_months[which(all_months < m)]
  prev_months <- tail(prev_months, 6)
  idx_j <- which(months %in% prev_months)
  if(length(idx_j) < 1) next
  
  fasta_i <- fasta[idx_i]
  fasta_j <- fasta[idx_j]
  
  # Combine for distance calculation
  combined <- c(fasta_i, fasta_j)
  alignment <- as.DNAbin(msaConvert(msa(combined, method="Muscle"), type="seqinr::alignment"))
  
  # Calculate pairwise distances
  dist_matrix <- dist.dna(alignment, model="N", pairwise.deletion=TRUE)
  
  # Extract distances: rows = i, cols = j
  n_i <- length(idx_i)
  n_j <- length(idx_j)
  dist_sub <- as.matrix(dist_matrix)[1:n_i, (n_i+1):(n_i+n_j)]
  
  mean_dist <- mean(dist_sub)
  median_dist <- median(dist_sub)
  
  antigenic_drift_results <- rbind(antigenic_drift_results,
    data.frame(month = m, mean_distance = mean_dist, median_distance = median_dist))
  
  cat(m, ": mean antigenic drift distance =", mean_dist, " median =", median_dist, "\n")
}

write.csv(antigenic_drift_results, "monthly_antigenic_drift.csv", row.names=FALSE)

# --- Detect significant antigenic drift events ---
drift_mean <- antigenic_drift_results$mean_distance
drift_z <- scale(drift_mean)
threshold <- 2  # Z-score > 2 视为显著漂移

important_events <- antigenic_drift_results[abs(drift_z) > threshold, ]
if(nrow(important_events) > 0) {
  cat("Significant antigenic drift detected at:\n")
  print(important_events)
  write.csv(important_events, "significant_antigenic_drift_events.csv", row.names=FALSE)
} else {
  cat("No significant antigenic drift events detected.\n")
}





# 添加安全的多序列比对函数
safe_msa <- function(sequences, max_seq = 10) {  # 增加到100条序列
  # 限制序列数量避免崩溃
  if(length(sequences) > max_seq) {
    cat("序列数量过多 (", length(sequences), "), 随机采样", max_seq, "条序列\n")
    sequences <- sequences[sample(length(sequences), max_seq)]
  } 
  res_align <- msa(sequences, method="ClustalOmega")
  return(res_align)
}

div_list <- foreach(m = all_months, .combine=rbind, .packages=c("msa", "ape", "pegas", "Biostrings")) %dopar% {
  idx <- which(months == m)
  if(length(idx) < 2) return(NULL) # 少于2个skip
  
  cat("处理月份:", m, ", 序列数量:", length(idx), "\n")
  
  fasta_sub <- fasta[idx]
  res_align <- safe_msa(fasta_sub)
  
  if(is.null(res_align)) {
    cat("月份", m, "比对失败，跳过\n")
    return(NULL)
  }
  
  tryCatch({
    alignment <- as.DNAbin(msaConvert(res_align, type="seqinr::alignment"))
    seg_sites <- length(seg.sites(alignment))
    haps <- haplotype(alignment)
    hd <- hap.div(haps)
    pi <- nuc.div(alignment)
    data.frame(
      month = m,
      n_seq = length(idx),
      seg_sites = seg_sites,
      hap_div = hd,
      nuc_div = pi,
      stringsAsFactors=FALSE
    )
  }, error = function(e) {
    cat("月份", m, "分析失败:", e$message, "\n")
    return(NULL)
  })
}
#write.csv(div_list, "monthly_diversity_stats_test.csv", row.names=FALSE)
# 移除NULL结果
div_list <- div_list[!is.null(div_list),]
write.csv(div_list, "monthly_diversity_stats.csv", row.names=FALSE)

