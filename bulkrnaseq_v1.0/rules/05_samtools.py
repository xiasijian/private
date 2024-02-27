########################### samtools_to_bam #################################
rule samtools_to_bam:
    input:
        sam_file=outdir+"/04_hisat2_mapped/output/"+"{sample}"+"_filt"+".sam"
    output:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    threads:
        2
    log:
        log_file=outdir+"/05_samtools/log/"+"{sample}"+"_filt"+"_samtools_to_bam.log"
    params:
        shell_script=outdir+"/05_samtools/log/"+"{sample}"+"filt"+"_samtools_to_bam.sh",
        software=config["software"]["samtools"],
        threads=config["samtools_para"]["thread"],
	conda_env=config["conda_env"]["samtools_env"]
    shell:
        """
        echo -e 'source activate {params.conda_env} \n{params.software} sort -o {output.bam_sorted_file} {input.sam_file} \n{params.software} index {output.bam_sorted_file}' >& {params.shell_script}
        
        sh {params.shell_script} > {log.log_file}
        """
########################### mapping stat #################################

rule samtools_stat:
    input:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    output:
        bam_stat=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted_mapped_stat.txt"
    params:
        shell_script=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted_mapped_stat.sh",
        software=config["software"]["samtools"]
    shell:
        """
        echo -e '{params.software} flagstat {input.bam_sorted_file} > {output.bam_stat}' >& {params.shell_script}
        
        {params.software} flagstat {input.bam_sorted_file} > {output.bam_stat}
        
        """
        
