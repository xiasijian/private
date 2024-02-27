"""
Author:
Sijian Xia

Date:
2023.09.06

reference paper:
(1) From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis, 2020
(2) Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data, 2013
(3) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia, 2012
(4) Methods for ChIP-seq analysis: A practical workflow and advanced applications, 2021

reference material:
(1) https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-5-chipseq/Epigenetics.html
(2) http://www.begenomics.com/
usage:
This snakemake pipeline is designed for chipseq or atac-seq

"""


import os
import sys
import glob
import pandas as pd
import re

#######################################################
## parameter setting
#######################################################

dir = str(os.getcwd())

outdir = dir + "/output/"
## create the output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

print("The workdirectory is ", dir)
configfile: os.path.join(dir,'config/config.yaml')


DATA_TYPE = config["run_para"]["data_type"]
DATA_DIR = config["run_para"]["data_dir"]
GENOME_FASTA_FILE  = config["run_para"]["genome"]
GFF_FILE = config["run_para"]["annotation_file"]

#######################################################
## Useful function
#######################################################

def create_condition_experiment_and_sample_dict(df,condition_select):
    df_filt = df[df['condition'].isin([condition_select])]
    sample_experiment = dict(zip(df_filt['sample'], df_filt['experiment']))
    return(sample_experiment)


def create_experiment_txt_file(experiment):
    text_file=open(outdir + "/"+experiment + ".txt","w")
    text_file.write(experiment)
    text_file.close()


#######################################################
## sample information
#######################################################
units  = pd.read_csv("sample_information.csv")
SAMPLES = units["sample"].unique().tolist()
EXPERIMENT = units["experiment"].unique().tolist()

get_fq1_from_experiment = units.set_index("experiment")['fq1'].to_dict()
get_fq2_from_experiment = units.set_index("experiment")['fq2'].to_dict()

get_treat_from_sample = create_condition_experiment_and_sample_dict(df=units,condition_select="IP")
get_control_from_sample = create_condition_experiment_and_sample_dict(df=units,condition_select="Input_control")


for i in units["experiment"]:
    print(i)
    if os.path.exists(outdir + "/"+i + ".txt"):
        print("files exited")
    else:
        create_experiment_txt_file(i)


#######################################################
## OUTPUT CHECK
#######################################################

rule all:
    input:
        r1_filt = expand(outdir + "/01_clean_fastq/" + "{experiment}" + "_val_1.fq.gz",experiment=EXPERIMENT),
        r2_filt = expand(outdir + "/01_clean_fastq/" + "{experiment}" + "_val_2.fq.gz",experiment=EXPERIMENT)
        

#######################################################
# Rules
#######################################################

include : "rules/01_clean_fastq.py"