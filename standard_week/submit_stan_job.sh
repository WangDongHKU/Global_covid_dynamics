#!/bin/bash

# SLURM提交脚本 - Stan模型运行
# 使用方法: sbatch submit_stan_job.sh

#SBATCH --job-name=stan_covid_model
#SBATCH --output=stan_job_%j.out
#SBATCH --error=stan_job_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com

# 设置工作目录
cd /home/dongw21/Global_COVID/standard_week

# 创建进度文件
PROGRESS_FILE="/home/dongw21/Global_COVID/standard_week/progress_log.txt"
echo "=== SLURM任务开始 ===" > "$PROGRESS_FILE"
echo "任务ID: $SLURM_JOB_ID" >> "$PROGRESS_FILE"
echo "开始时间: $(date)" >> "$PROGRESS_FILE"
echo "节点: $SLURM_NODELIST" >> "$PROGRESS_FILE"
echo "CPU核心数: $SLURM_CPUS_PER_TASK" >> "$PROGRESS_FILE"
echo "" >> "$PROGRESS_FILE"

# 运行R脚本
echo "开始运行Stan模型..." >> "$PROGRESS_FILE"
Rscript standard_model_10.R

# 记录任务完成
echo "" >> "$PROGRESS_FILE"
echo "=== SLURM任务完成 ===" >> "$PROGRESS_FILE"
echo "完成时间: $(date)" >> "$PROGRESS_FILE"
echo "任务ID: $SLURM_JOB_ID" >> "$PROGRESS_FILE"

echo "任务完成！进度文件: $PROGRESS_FILE"
