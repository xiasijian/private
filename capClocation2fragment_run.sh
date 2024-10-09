#!/bin/bash

source activate /lustre/home/sjxia/miniconda3/install/envs/capc_map

restfragfile=/lustre/home/sjxia/04_ref_by_species/ucsc/mm9/capture_c/Dpnii/dpnII_map_mm9.bed
outfile=target_probe24_dpnii.bed
inputfile=target_probe24.bed
capClocation2fragment -r ${restfragfile} -o ${outfile} -i ${inputfile}
