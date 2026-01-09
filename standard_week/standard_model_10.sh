#!/bin/bash
#SBATCH --job-name=COVID_S        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL   # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --ntasks=4                   # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --mem=130G                     # 6. Request total amount of RAM
#SBATCH --time=15-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=/scr/u/dongw21/Global_COVID/standard_week/%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=/scr/u/dongw21/Global_COVID/standard_week/%x_%j.err             #    Standard error log as $job_name_$job_id.err

# 创建进度文件
module load cmake
module load gcc
module load R/4.4.3
mkdir -p /home/dongw21/R_libs/4.4.3
Rscript /scr/u/dongw21/Global_COVID/standard_week/standard_model_10.R
