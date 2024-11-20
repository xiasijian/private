#!/bin/bash
#SBATCH --job-name=hic_run

#SBATCH --partition=cpuPartition

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=8

#SBATCH --error=%j.err
#SBATCH --output=%j.out

module load anaconda
source activate cooltools_env

cooler merge GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.cool  GSE247254_YY1aid_async_untreated.mcool::resolutions/5000 GSE254373_LDB1_AID_MICRO_C_untreated_merged.mcool::resolutions/5000

module load anaconda
source activate hictk

cooler zoomify --balance -r 5000 -n 8 GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.cool -o GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.cool
hictk convert -t 8 GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.cool GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.hic

cooltools expected-cis -p 16 -o GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.expected_cis.tsv --view mm9_chr1toX.bed --clr-weight-name weight GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.cool::/resolutions/5000


cooltools dots --view mm9_chr1toX.bed -p 10 --max-loci-separation 2000000 --fdr 0.02 --clustering-radius 39000 -o GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.dot.txt GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.cool::/resolutions/5000 GSE254373_Ldb1_GSE247254_Yy1_aid_microc_async_untreated_merge.5000.balance.expected_cis.tsv

