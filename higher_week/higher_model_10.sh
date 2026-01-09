#!/bin/bash
#SBATCH --job-name=sCOVI-H        # 1. Job name
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # 申请 4 个 CPU 线程
#SBATCH --mem=100G                  # 内存（根据情况调）
#SBATCH --time=6-10:00:00            # 最长运行时间 6 天 10 小时
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --output=/scr/u/dongw21/Global_COVID/higher_week/%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=/scr/u/dongw21/Global_COVID/higher_week/%x_%j.err             #    Standard error log as $job_name_$job_id.err

# 加载必要的模块
module load cmake
module load gcc
module load root
module load R/4.4.3
mkdir -p /home/dongw21/R_libs/4.4.3
Rscript //scr/u/dongw21/Global_COVID/higher_week/Higher_model_10.R
