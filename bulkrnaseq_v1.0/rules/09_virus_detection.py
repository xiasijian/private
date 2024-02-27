rule extract_unmapped:
    input:
        bam_sorted_file=outdir+"/05_samtools/output/"+"{sample}"+"_filt"+"_sorted.bam"
    output:
        unmapped_file=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped.fastq"
    log:
        log_file=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_extract_umapped.sh"
    params:
        samtools=config["software"]["samtools"]
    shell:
        """
        echo -e '{params.samtools} fastq -f 4 {input.bam_sorted_file} > {output.unmapped_file}' > {log.log_file}
        
        {params.samtools} fastq -f 4 {input.bam_sorted_file} > {output.unmapped_file}
        
        """
        
rule bowtie2_mapped:
    input:
        unmapped_file=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped.fastq"
    output:
        umapped_bam=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped_bowtie2.bam"
    log:
        logfile=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped_bowtie2.log"
    params:
        bowtie2_env=config["conda_env"]["bowtie2_env"],
        umapped_sam=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped_bowtie2.sam",
        subject_database=config["virus_detection_para"]["subject_databse"],
        shell_script=outdir+"/09_virus_detection/output/"+"{sample}"+"_filt_sort_umapped_bowtie2.sh",
        samtools_env=config["conda_env"]["samtools_env"]
        
    shell:
        """
        echo -e 'source activate {params.bowtie2_env} \nbowtie2 -x {params.subject_database} -U {input.unmapped_file} -p 4 -S {params.umapped_sam} \nsource activate {params.samtools_env} \nsamtools view -@4 -b {params.umapped_sam} > {output.umapped_bam}' > {params.shell_script}
        
        sh {params.shell_script} >& {log.logfile}
        """
    