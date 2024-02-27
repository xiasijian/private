rule RSeQC_read_distribution:
    input:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    output:
        bam_qc_read_distribution=outdir+"/06_RSeQC/output/"+"{sample}"+"_filt"+"_sorted_read_distribution.txt"
    threads:
        2
    params:
        read_distribution_run=config["software"]["RSeQC"]+"/bin/read_distribution.py",
        bed_file=config["RSeQC_para"]["bed_file"],
        RSeQC_env=config["software"]["RSeQC"]
        
    shell:
        """
        source activate {params.RSeQC_env}
        {params.read_distribution_run} -i {input.bam_sorted_file} -r {params.bed_file} > {output.bam_qc_read_distribution}
    
    
        """