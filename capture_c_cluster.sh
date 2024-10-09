#!/bin/bash

#SBATCH --job-name=mcc_run4
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --error=%j.err
#SBATCH --output=%j.out


data_dir=/lustre/home/sjxia/03_project/3d_genome_seq_trainning/03_capture_c/experiment/capturec_mcc_20240926data/01_rawdata/arrange_data/

cat sample_list.txt | parallel --jobs 2 sh capture_c_pipeline.sh ${data_dir} {}
