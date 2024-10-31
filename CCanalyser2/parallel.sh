#!/bin/bash

#SBATCH --job-name=mcc_run
#SBATCH --partition=cpuPartition
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=10
#SBATCH --error=%j.err
#SBATCH --output=%j.out

data_dir=/lustre/home/sjxia/03_project/3d_genome_seq_trainning/03_capture_c/experiment/lyz_data20241012/02_analysis/v4
cat sample_list.txt | parallel -j 3 sh run_pipeline.sh {} ${data_dir}

