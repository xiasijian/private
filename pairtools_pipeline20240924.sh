#!/bin/bash

sample_name=$1
data_dir=$2

input_r1=${data_dir}/${sample_name}_raw_1.fq.gz
input_r2=${data_dir}/${sample_name}_raw_2.fq.gz

module load anaconda
source activate chipseq_env


## https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html#typical-hi-c-workflow

##------------------------------
## mapping ot mm9
##------------------------------

genome_prefix=mm9
bwa_index=/lustre/home/sjxia/04_ref_by_species/ucsc/${genome_prefix}/index/bwa/${genome_prefix}.fa

sam_file=${sample_name}_bwa_${genome_prefix}_mapped.sam
thread=7
# bwa mem -t ${thread} -v 3 -SP ${bwa_index} ${input_r1} ${input_r2} > ${sam_file}

##------------------------------
## sam to bam
##------------------------------

samtools=~/.conda/envs/samtools_env/bin/samtools
bam_sort=${sample_name}_bwa_${genome_prefix}_mapped_sorted.bam
tmp_bam=${sample_name}.tmp.bam
thread=7
# ${samtools} sort --threads ${thread} -T ${tmp_bam} ${sam_file} -o ${bam_sort}


##------------------------------
## pairtools parse
##------------------------------
## finds ligation junctions and records the outer-most 5' aligned base pair and 
## the strand of each one of the paired reads into a .pairsam file

source activate hic_env

chroms_file=/lustre/home/sjxia/04_ref_by_species/ucsc/${genome_prefix}/genome/${genome_prefix}.genome
pairsam=${sample_name}_bwa_${genome_prefix}_mapped.pairsam


# pairtools parse --walks-policy 5unique --add-columns mapq -c ${chroms_file} --assembly ${genome_prefix} ${sam_file} > ${pairsam}

##------------------------------
## pairtools sort
##------------------------------
tmp_dir=./
sort_pairsam=${sample_name}_bwa_${genome_prefix}_mapped_sorted.pairsam
# pairtools sort --nproc 5 --tmpdir=${tmp_dir} ${pairsam} > ${sort_pairsam}


dedup_pairsam=${sample_name}_bwa_${genome_prefix}_mapped_sorted_dedup.pairsam
# pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats ${sample_name}_pcr_duplicated_stats.txt --output ${dedup_pairsam} ${sort_pairsam}

##------------------------------
## pairtools split
##------------------------------
mapped_dedup_pairs=${sample_name}_bwa_${genome_prefix}_mapped_sorted_dedup.pairs
split_bam=${sample_name}_bwa_${genome_prefix}_mapped_sorted_dedup_split.bam

# pairtools split --nproc-in 8 --nproc-out 8 --output-pairs ${mapped_dedup_pairs} --output-sam ${split_bam} ${dedup_pairsam}



##------------------------------
## convert pair to hic
##------------------------------
mapped_dedup_hic_file=${sample_name}_bwa_${genome_prefix}_mapped_sorted_dedup.hic
pair_to_hic_log=${sample_name}_bwa_${genome_prefix}_mapped_sorted_dedup.log
java -Xmx5g -jar /lustre/home/sjxia/02_software/juicer_tools.jar pre ${mapped_dedup_pairs} ${mapped_dedup_hic_file} ${genome_prefix} >& ${pair_to_hic_log}


##------------------------------
## Load pairs to cooler
##------------------------------
# source activate hic_env

# bgzip -@ 5 -c ${mapped_dedup_pairs} > ${mapped_dedup_pairs}.gz
# pairix ${mapped_dedup_pairs}.gz


## Binning and matrix formats
# Let’s convert our .pairs file to a .cool file binned at 10Kb resolution. When binning the matrix a high resolution bin size (1kb - 10kb) is recommended because a lower resolution can be easily obtained by summing adjacent bins.
#cooler cload pairix -p 16 ${chroms_file}:1000 ${mapped_dedup_pairs}.gz ${sample_name}.matrix_1kb.cool


## Genereting multi-resolutions files and visualizing the contact matrix
# cooler zoomify --balance -p 16 ${sample_name}.matrix_1kb.cool