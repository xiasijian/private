rule high_mapped:
    input:
        sorted_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam",
        sorted_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_sorted.bam.bai"
    output:
        mapped_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30.bam",
        mapped_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30.bam.bai"
    log:
        outdir + "/02_align/" + "{experiment}"+".mapped.q30.bam.log"
    params:
        samtools = config["software"]["samtools"],
        threads = config["samtools_para"]["threads"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_mapped_q30_bam.sh"
    shell:
        """
        rm -f {params.shell_script}
        echo -e '{params.samtools} view -bF 1804 -q 30 -@ {params.threads} {input.sorted_bam} > {output.mapped_bam}' > {params.shell_script}
        echo -e '{params.samtools} index {output.mapped_bam}' >> {params.shell_script}
        bash {params.shell_script} > {log}
        
        """


rule rmdp:
    input:
        mapped_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30.bam",
        mapped_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30.bam.bai"
    output:
        rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam",
        rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam.bai"
    log:
         outdir + "/02_align/" + "{experiment}"+".sorted.rmdup.bam.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        picard = config["software"]["picard"],
        samtools = config["software"]["samtools"],
        pca_rm = config["picard_para"]["pca_rm"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_picard_rmdp.sh",
        metrics_file = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_metrics.txt"
    threads: 2 
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'source activate {params.chipseq_env}' >> {params.shell_script}
        echo -e '{params.picard} MarkDuplicates -I {input.mapped_bam} {params.pca_rm} -O {output.rmdp_bam} --METRICS_FILE {params.metrics_file}' >> {params.shell_script}
        echo -e '{params.samtools} index {output.rmdp_bam}' >> {params.shell_script}
        bash {params.shell_script} > {log}


        """


#  -v      Only report those entries in A that have _no overlaps_ with B.
#             - Similar to "grep -v" (an homage).

rule rmblacklist:
    input:
        rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam",
        rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam.bai"
    output:
        blacklist_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam",
        blacklist_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam.bai"
    log:
        outdir + "/02_align/" + "{experiment}"+".sorted.rmdup.rmblacklist.bam.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        samtools = config["software"]["samtools"],
        blacklist = config["run_para"]["blacklist"],
        shell_script = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_picard_rmdp_rmblacklist.sh"
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'module load anaconda \nsource activate {params.chipseq_env}' >> {params.shell_script}
        echo -e 'intersectBed -v -abam {input.rmdp_bam} -b {params.blacklist} > {output.blacklist_bam}' >> {params.shell_script}
        echo -e '{params.samtools} index {output.blacklist_bam}' >> {params.shell_script}
        bash {params.shell_script} > {log}
        """