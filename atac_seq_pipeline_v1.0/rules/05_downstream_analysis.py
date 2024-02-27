rule deeptools_bamcoverage_set1:
    input:
        blacklist_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam",
        blacklist_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam.bai"
    output:
        cpm_bw = outdir + "/04_downstream/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist_cpm.bw"
    log:
        outdir + "/04_downstream/" + "{experiment}" + "_deeptools_bamcoverage_cpm_set1.log"
    params:
        deeptools_env = config["conda_env"]["deeptools_env"],
        genomesize = str(config['macs2_para']['genomesize']),
        set1 = str(config['deeptools']['bamcoverage']["set1"]),
        shell_script = outdir + "/04_downstream/" + "{experiment}" + "_deeptools_bamcoverage_cpm_set1.sh"
    threads:
        4
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'module load anaconda \nsource activate {params.deeptools_env}' > {params.shell_script}
        echo -e 'bamCoverage -of bigwig {params.set1} -b {input.blacklist_bam} \\
            -o {output.cpm_bw} -p {threads}' >> {params.shell_script}
        
        bash {params.shell_script}
        
        """

