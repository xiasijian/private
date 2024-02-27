## 官方教程
## https://subread.sourceforge.net/featureCounts.html
## https://subread.sourceforge.net/RNAseqCaseStudy.html
## 参数说明
#解释：-T 1：线程数为1；-p 表示双端测序 -t exon表示 feature名称为exon；-g gene_id表示meta-feature名称为gene_id(ensembl名称)；-a 表示输入的GTF基因组注释文件；-o 设置输出文件类型
# 1、feature是指基因组上被定义的一个片段区域，meta-feature是指多个feature组成的区域，如exon和gene的关系 2、分享相同的feature identifier（GTF文件中有） 的features属于同一个meta-feature 3、featurecounts可以对features和meta-feature进行计数
# -f：该参数设置后统计的是feature层面（默认就是exon）的参数，如果不设置则是直接统计meta-feature参数（及一个gene中的多个exon）
# -f If specified, read summarization will be performed at feature
# level (eg. exon level). Otherwise, it is performed at metafeature level (eg. gene level).

########################### featurecount #################################
rule featureCounts_gene:
    input:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    output:
        gene_count_file=outdir+"/07_featureCounts_gene/output/"+"{sample}"+"_featureCounts_gene_f_"+"output.txt"
    threads:
        4
    log:
        log_file=outdir+"/07_featureCounts_gene/log/"+"{sample}"+"_filt"+"_featureCounts_gene_f.log"
    params:
        shell_script=outdir+"/07_featureCounts_gene/log/"+"{sample}"+"_filt"+"_featureCounts_gene_f.sh",
        software=config["software"]["featureCounts"],
        ref_file=config["featureCounts_para"]["ref_file"],
        meta_feature=config["featureCounts_para"]["meta_feature"],
        feature_name=config["featureCounts_para"]["feature_name"],
        threads=config["featureCounts_para"]["threads"]
    shell:
        """
        echo -e '{params.software} -p -T {params.threads} -t gene -f -g {params.meta_feature} -o {output.gene_count_file} -a {params.ref_file} {input.bam_sorted_file}' > {params.shell_script}
        
        {params.software} -p -T {params.threads} -t gene -f -g {params.meta_feature} -o {output.gene_count_file} -a {params.ref_file} {input.bam_sorted_file} >& {log.log_file}
        
        """
        
