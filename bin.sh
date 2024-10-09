#!/bin/bash

source activate /lustre/home/sjxia/miniconda3/install/envs/capc_map

sample_name=cap-c-FA-max-G4-0h-iaa
species=mm9
gene=Max_C1
save_name=${sample_name}_capturec_probe24
chrom_size=/lustre/home/sjxia/04_ref_by_species/ucsc/${species}/genome/${species}.genome

cis_raw_bdg=${save_name}_rawpileup_${gene}.bdg
## out cis
bin=10
windows=1000
cis_norm_bdg=${save_name}_bin${bin}_win${windows}_${gene}.bdg
total_interaction_from_report=192802


capCpileup2binned -i ${cis_raw_bdg} -o ${cis_norm_bdg} -c ${chrom_size} -t ${sample_name}_${gene} -b ${bin} ${windows} -n ${total_interaction_from_report}

sample_name=cap-c-FA-max-G4-24h-iaa
species=mm9
gene=Max_C1
save_name=${sample_name}_capturec_probe24
chrom_size=/lustre/home/sjxia/04_ref_by_species/ucsc/${species}/genome/${species}.genome

cis_raw_bdg=${save_name}_rawpileup_${gene}.bdg
## out cis
bin=10
windows=1000
cis_norm_bdg=${save_name}_bin${bin}_win${windows}_${gene}.bdg
total_interaction_from_report=218944


capCpileup2binned -i ${cis_raw_bdg} -o ${cis_norm_bdg} -c ${chrom_size} -t ${sample_name}_${gene} -b ${bin} ${windows} -n ${total_interaction_from_report}

