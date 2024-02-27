
import pandas as pd
import numpy as np
import os
import sys
import glob

## get the current directory
dir=str(os.getcwd())
configfile: os.path.join(dir,'config/config.yaml')

##########################################
## Sample information
##########################################
sample_file=dir+"/"+"sample_info.csv"
samples=pd.read_csv(sample_file,sep=",",keep_default_na=False, na_values=[""])
sample_list=list(samples["sample"].unique())
sample_list = [x for x in sample_list if str(x) != 'nan']

if not os.path.exists(dir+"/output/"):
    os.makedirs(dir+"/output/") 

outdir=dir+"/output"
data_dir=config["run_para"]["data_dir"]
fastq_suffix=config["run_para"]["fastq_suffix"]

rule all:
    input:
        ## fastqc
        #fastqc_raw=expand(outdir+"/01_fastqc_raw/"+"output/"+"{sample}"+"fastq_raw_completed.txt",sample=sample_list),
        #fastqc_filt=expand(outdir+"/03_fastqc_filt/"+"output/"+"{sample}"+"fastq_filt"+"_completed.txt",sample=sample_list),
        ## fastp
        filt_R1=expand(outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R1_filt"+".fq.gz",sample=sample_list),
        filt_R2=expand(outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R2_filt"+".fq.gz",sample=sample_list),
        json=expand(outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_filt"+"_fastp.json",sample=sample_list),
        ## hisat2
        #sam_file=expand(outdir+"/04_hisat2_mapped/output/"+"{sample}"+"_filt"+".sam",sample=sample_list),
        ## samtools
        bam_sorted_file=expand(outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam",sample=sample_list),
        bam_stat=expand(outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted_mapped_stat.txt",sample=sample_list),
        ## RSeQC
	#bam_qc_read_distribution=expand(outdir+"/06_RSeQC/output/"+"{sample}"+"_filt"+"_sorted_read_distribution.txt",sample=sample_list),
        ## featureCounts
        #count_file=expand(outdir+"/07_featureCounts/output/"+"{sample}"+"_featureCounts_"+"output.txt",sample=sample_list),
        ## featureCounts change parameter
        gene_count_file=expand(outdir+"/07_featureCounts_gene/output/"+"{sample}"+"_featureCounts_gene_f_"+"output.txt",sample=sample_list),
        ## htseq-count
        htseq_count_file=expand(outdir+"/08_htseq_count/output/"+"{sample}"+"_htseq_count_gene_"+"output.txt",sample=sample_list),
        ## virus detection
        #umapped_bam=expand(outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped_bowtie2.bam",sample=sample_list)

        
include: "rules/02_fastp_filt.py"

include: "rules/04_hisat2_mapped.py"
include: "rules/05_samtools.py"
#include: "rules/06_RSeQC.py"
#include: "rules/07_featureCount.py"

## featureCount change parameter
include: "rules/07_featureCount_gene.py"

## htseq-count
include: "rules/08_htseq-count.py"

## virus detection
#include: "rules/09_virus_detection.py"

#include: "rules/01_fastqc_raw_qc.py"
#include: "rules/03_fastqc_filt_qc.py"

