rule star_human:
    input:
        filt_R1=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R1_filt"+".fq.gz",
        filt_R2=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R2_filt"+".fq.gz"
    output:
        unmap_R1=outdir+"/11_venus_stimulate/"+"output/"+"{sample}"+"/"+"Unmapped.out.mate1",
        unmap_R2=outdir+"/11_venus_stimulate/"+"output/"+"{sample}"+"/"+"Unmapped.out.mate2"
    log:
        log_file=outdir+"/11_venus_stimulate/"+"{sample}"+"_run_star_mapped.log"
    threads: 5
    params:
        venus_env=config["conda_env"]["venus_env"],
        shell_script=outdir+"/11_venus_stimulate/"+"{sample}"+"_run_star_mapped.sh",
        humanGenome="/share/home/xiasj/01_biosoft/Venus/indices/human2.genomeDir",
        outdir=outdir+"/11_venus_stimulate/output/"+"{sample}"+"/"
    shell:
        """
        echo -e 'source activate {params.venus_env} \nSTAR --runThreadN 5 --outFileNamePrefix {params.outdir} --genomeDir {params.humanGenome} --readFilesCommand zcat --readFilesIn {input.filt_R1} {input.filt_R2} --outSAMtype None --outReadsUnmapped Fastx' > {params.shell_script}
        
        sh {params.shell_script} >& {log.log_file}
        
        """