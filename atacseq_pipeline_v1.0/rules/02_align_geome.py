rule align:
    input:
        r1_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_1.fq.gz",
        r2_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_2.fq.gz",
        clean_fastq_completed = outdir + "/01_clean_fastq/output/" + "{experiment}" +"_completed.txt"
        #r1_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "/"+"{experiment}" +"_1.fastq.gz",
        #r2_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "/"+"{experiment}" +"_2.fastq.gz",
    output:
        sam_file =  outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped.sam"
    log:
        outdir + "/02_align/" + "{experiment}" + "_bowtie2_align.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        bowtie2 = config["software"]["bowtie2"],
        #bowtie2_para = " ".join(config["bowtie2_para"].values()), #take argument separated as a list separated with a space
        bowtie2_para = "-X2000 --no-mixed --no-discordant -k1 -N1 " + config["bowtie2_para"]["index"],
        threads = config["samtools_para"]["threads"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_bowtie2_align.sh"
    threads: config["samtools_para"]["threads"]

    shell:
        """
        rm -f {params.shell_script}
        echo -e 'source activate {params.chipseq_env}' >> {params.shell_script}
        echo -e '{params.bowtie2} {params.bowtie2_para} --threads {params.threads} -1 {input.r1_filt} -2 {input.r2_filt} -S {output.sam_file}' >> {params.shell_script}

        bash {params.shell_script} >& {log}

        """
rule sam_to_bam:
    input:
        sam_file =  outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped.sam"
    output:
        sorted_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam",
        sorted_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam.bai"
    log:
        outdir + "/02_align/" + "{experiment}" + "_samtool_sam_to_bam.log"
    params:
        samtools = config["software"]["samtools"],
        threads = config["samtools_para"]["threads"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_samtool_sam_to_bam.sh"
    threads: config["samtools_para"]["threads"]
    shell:
        """
        rm -f {params.shell_script}
        echo -e '{params.samtools} sort -@ {params.threads} -O bam -o {output.sorted_bam} {input.sam_file}' >> {params.shell_script}
        echo -e 'rm -f {input.sam_file}' >> {params.shell_script}
        echo -e '{params.samtools} index {output.sorted_bam}' >> {params.shell_script}

        bash {params.shell_script} >& {log}

        """