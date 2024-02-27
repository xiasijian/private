## 官方教程
## https://subread.sourceforge.net/featureCounts.html

## 这里可能有一个错误，就是没有指定-p，表示数据为paired-end，双末端测序数据
## 参数说明
#解释：-T 1：线程数为1；-t exon表示 feature名称为exon；-g gene_id表示meta-feature名称为gene_id(ensembl名称)；-a 表示输入的GTF基因组注释文件；-o 设置输出文件类型

########################### featurecount #################################
rule featureCounts:
    input:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    output:
        count_file=outdir+"/07_featureCounts/output/"+"{sample}"+"_featureCounts_"+"output.txt"
    threads:
        2
    log:
        log_file=outdir+"/07_featureCounts/log/"+"{sample}"+"_filt"+"_featureCounts.log"
    params:
        shell_script=outdir+"/07_featureCounts/log/"+"{sample}"+"_filt"+"_featureCounts.sh",
        software=config["software"]["featureCounts"],
        ref_file=config["featureCounts_para"]["ref_file"],
        meta_feature=config["featureCounts_para"]["meta_feature"],
        feature_name=config["featureCounts_para"]["feature_name"],
        threads=config["featureCounts_para"]["threads"]
    shell:
        """
        echo -e '{params.software} -T {params.threads} -t {params.feature_name} -g {params.meta_feature} -o {output.count_file} -a {params.ref_file} {input.bam_sorted_file}' > {params.shell_script}
        
        {params.software} -T {params.threads} -t {params.feature_name} -g {params.meta_feature} -o {output.count_file} -a {params.ref_file} {input.bam_sorted_file} >& {log.log_file}
        
        """
        
