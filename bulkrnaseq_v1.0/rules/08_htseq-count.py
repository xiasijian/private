## 官方文档
## https://htseq.readthedocs.io/en/master/htseqcount.html
## 参数说明：-r pos 表示按照比对位置进行排序; -n是指定CPU数；-a Skip all reads with MAPQ alignment quality lower than the given minimum value (default: 10). -f表示输入的文件格式

rule htseq_count:
    input:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    output:
        htseq_count_file=outdir+"/08_htseq_count/output/"+"{sample}"+"_htseq_count_gene_"+"output.txt"
    threads: 4
    log:
        log_file=outdir+"/08_htseq_count/log/"+"{sample}"+"_filt_"+"htseq_count_gene.log"
    params:
        shell_script=outdir+"/08_htseq_count/log/"+"{sample}"+"_filt_"+"htseq_count_gene.sh",
        software=config["software"]["htseq_count"],
        ref_file=config["featureCounts_para"]["ref_file"]
    shell:
        """
        echo -e '{params.software} -r pos -n 4 -a 10 -f bam {input.bam_sorted_file} {params.ref_file} > {output.htseq_count_file}' > {params.shell_script}
        
        {params.software} -r pos -n 4 -a 10 -f bam {input.bam_sorted_file} {params.ref_file} > {output.htseq_count_file} 
        
        """