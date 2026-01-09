#!/bin/bash
#SBATCH --job-name=Visualization        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL   # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --ntasks=1                   # 5. Request total number of tasks (MPI workers)
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --mem=80G                     # 6. Request total amount of RAM
#SBATCH --time=9-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=/home/dongw21/Global_COVID/standard_week/%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=/home/dongw21/Global_COVID/standard_week/%x_%j.err             #    Standard error log as $job_name_$job_id.err

# åå»ºè¿åº¦æä»¶
PROGRESS_FILE="/home/dongw21/Global_COVID/standard_week/progress_log.log"
module load cmake
module load gcc
module load R/4.4.3
#mkdir -p /home/dongw21/R_libs/4.4.3

Rscript /home/dongw21/Global_COVID/standard_week/data_analysis_visualision.R
