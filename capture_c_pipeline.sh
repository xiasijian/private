#!/bin/bash

data_dir=$1
sample_name=$2

module load anaconda
source activate mcc_env
##--------------------------------
## trim_galore
##--------------------------------
input_fq1=${data_dir}/${sample_name}_raw_1.fq.gz
input_fq2=${data_dir}/${sample_name}_raw_2.fq.gz

threads=4
outdir=./
trim_galore_soft=/lustre/home/sjxia/02_software/TrimGalore-0.6.10/trim_galore
#${trim_galore_soft} --cores ${threads} --trim-n --paired --output_dir ${outdir} ${input_fq1} ${input_fq2} >& ${sample_name}_trim.log


##--------------------------------
## digest fastq
##--------------------------------
source activate /lustre/home/sjxia/miniconda3/install/envs/capc_map

trim_fq1=${sample_name}_raw_1_val_1.fq
trim_fq2=${sample_name}_raw_2_val_2.fq


# gunzip -c ${trim_fq1}.gz > ${trim_fq1}
# gunzip -c ${trim_fq2}.gz > ${trim_fq2}

enzyme=dpnii
# digest_fastq=${sample_name}_trim_digested${enzyme}.fastq
# rm -f ${digest_fastq}

# capCdigestfastq -1 ${trim_fq1} -2 ${trim_fq2} -e GATC -p 5 -o ${digest_fastq}


##--------------------------------
## mapping to mm9
##--------------------------------
source activate /lustre/home/sjxia/miniconda3/install/envs/capc_map
species=mm9
ref_index=/lustre/home/sjxia/04_ref_by_species/ucsc/${species}/index/bowtie/${species}
align_sam=${sample_name}_trim_digested${enzyme}_bowtie_${species}.sam

# bowtie -p 5 -m 1 --best --strata --sam --chunkmb 256 ${ref_index} ${digest_fastq} ${align_sam}

##--------------------------------
## samtools process
##--------------------------------
bam_file=${sample_name}_trim_digested${enzyme}_bowtie_${species}.bam
bam_sort=${sample_name}_trim_digested${enzyme}_bowtie_${species}_sorted.bam
sam_sort=${sample_name}_trim_digested${enzyme}_bowtie_${species}_sorted.sam

# samtools view -S -b -@ 12 -o ${bam_file} ${align_sam}
# samtools sort -n -@ 12 -o ${bam_sort} -T tempsort ${bam_file}
# samtools view -h -@ 12 -o ${sam_sort} ${bam_sort}

genome_enzyme=/lustre/home/sjxia/04_ref_by_species/ucsc/mm9/capture_c/Dpnii/dpnII_map_mm9.bed
target_bed=target_probe24_dpnii.bed
## need to change
save_name=${sample_name}_capturec_probe24

#capCmain -r ${genome_enzyme} -t ${target_bed} -s ${sam_sort} -o ${save_name} -e 1000 -i


target_loc=chr12:78062681-78063197
gene=Max_C1
cis_pairs=${save_name}_validpairs_${gene}.pairs
cis_raw_bdg=${save_name}_rawpileup_${gene}.bdg

trans_pairs=${save_name}_validinterchom_${gene}.pairs
trans_raw_bdg=${save_name}_rawpileup_interchom_${gene}.bdg



# capCpair2bg -i ${cis_pairs} -o ${cis_raw_bdg} -n ${gene} -t ${target_loc}
# capCpair2bg -i ${trans_pairs} -o ${trans_raw_bdg} -n ${gene} -t ${target_loc} --interchrom


