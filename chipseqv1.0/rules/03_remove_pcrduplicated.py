rule rmdp:
    input:
        sorted_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam",
        sorted_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam.bai"
    output:
        rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted_rmdp.bam",
        rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted_rmdp.bam.bai"
    log:
         outdir + "/02_align/" + "{experiment}"+".sorted.rmdup.bam.log"
    threads: 2 
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        picard = config["software"]["picard"],
        samtools = config["software"]["samtools"],
        pca_rm = config["picard_para"]["pca_rm"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted_picard_rmdp.sh"
        metrics_file = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted_rmdp_metrics.txt"
    shell:
    """
    rm -f {params.shell_script}
    echo -e 'module load anaconda \nsource activate {params.chipseq_env}' >> {params.shell_script}
    echo -e '{params.picard} MarkDuplicates -I {input.sorted_bam} {params.pca_rm} -O {params.rmdp_bam} --METRICS_FILE {params.metrics_file}' >> {params.shell_script}
    echo -e '{params.samtools} index {output.rmdp_bam}' >> {params.shell_script}
    sh {params.shell_script} > {log}
    
     """
        