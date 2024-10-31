#!/bin/bash

#SBATCH --job-name=mcc_run
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --error=%j.err
#SBATCH --output=%j.out


module load anaconda
source activate mcc_env

sample_name=$1
data_dir=$2

## data_dir=/lustre/home/sjxia/03_project/3d_genome_seq_trainning/03_capture_c/experiment/lyz_data20241012/02_analysis/v4
## /lustre/home/sjxia/03_project/3d_genome_seq_trainning/03_capture_c/CCanalyser2/dpngenome3_1.pl ../../genome/mm9.fa

##---------------------------------------------
## flash
##---------------------------------------------

flash --interleaved-output -t 10 -M 150 ${sample_name}_raw_1_val_1.fq ${sample_name}_raw_2_val_2.fq --output-prefix ${sample_name}_flash >& ${sample_name}_flash.log

flash_fq=${sample_name}_trim_flash.fastq

cat ${sample_name}_flash.notCombined.fastq ${sample_name}_flash.extendedFrags.fastq > ${flash_fq}


# ## Use of uninitialized value $hash{"new name"} in concatenation (.) or string at ./dpnII2E.pl line 67, <FH> line 31392.

##---------------------------------------------
## dpnII2E
##---------------------------------------------
./dpnII2E_v2.pl ${flash_fq}


index=/lustre/home/sjxia/04_ref_by_species/ucsc/mm9/index/bowtie/mm9


##---------------------------------------------
## bowtie
##---------------------------------------------
enzyme_cut_fastq=${sample_name}_trim_flash_REdig.fastq
sam=${data_dir}/${sample_name}_trim_flash_REdig.sam

~/miniconda3/install/envs/capc_map/bin/bowtie -p 10 -m 2 --best --strata --sam --chunkmb 256 -x ${index} -q ${enzyme_cut_fastq} --sam ${sam}


##---------------------------------------------
## CCanalyser3_v2
##---------------------------------------------

## mm9_sizes.txt
chrom_size_dir=/lustre/home/sjxia/04_ref_by_species/ucsc/mm9/genome
oligo_file=${data_dir}/Max_DpnII_CCanalyser.txt
outdir=./
genome_rz=/lustre/home/sjxia/04_ref_by_species/ucsc/mm9/capture_c/Dpnii_ccanalyser2/mm9_dpnII_coordinates.txt
genome_name=mm9


perl CCanalyser3_v2.pl -b ${chrom_size_dir} -f ${sam} -r ${genome_rz} --genome ${genome_name} \
	-o ${oligo_file} -s ${sample_name}_ccanalyser --pf ${outdir}
