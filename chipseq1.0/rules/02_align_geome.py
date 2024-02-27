rule align:
    input:
        r1_filt = outdir + "/01_clean_fastq/" + "{experiment}" + "_val_1.fq.gz",
        r2_filt = outdir + "/01_clean_fastq/" + "{experiment}" + "_val_2.fq.gz"
    output:
        sorted_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam",
        sorted_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam.bai"
    log:
        outdir + "/02_align/" + "{experiment}" + "_bowtie2_align.log"
        
    threads: config["samtools_para"]["threads"]
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        bowtie2 = config["software"]["bowtie2"],
        bowtie2_para = " ".join(config["bowtie2_para"].values()), #take argument separated as a list separated with a space
        samtools = config["software"]["samtools"],
        sam_file = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped.sam",
        threads = config["samtools_para"]["threads"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_bowtie2_align.sh"
    shell:
    """
    rm -f {params.shell_script}
    echo -e 'module load anaconda \nsource activate '
    echo -e '{params.bowtie2} {params.bowtie2_para} --threads {params.threads} -1 {input.r1_filt} -2 {input.r2_filt} -S {params.sam_file} ' >> {param.shell_script}
    echo -e '{params.samtools} sort -@ {params.threads} -O bam -o {output.sorted_bam} {params.sam_file}' >> {param.shell_script}
    echo -e 'rm -f {params.sam_file}' >> {param.shell_script}
    echo -e '{params.samtools} index {output.sorted_bam}' >> {params.shell_script}
    sh {params.shell_script} >& {log}
    
    """