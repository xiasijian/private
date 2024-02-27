########################### hista2 mapping #################################
rule hisat2_mapped:
    input:
        filt_R1=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R1_filt"+".fq.gz",
        filt_R2=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R2_filt"+".fq.gz",
        json=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_filt"+"_fastp.json"
    output:
        sam_file=outdir+"/04_hisat2_mapped/output/"+"{sample}"+"_filt"+".sam"
    threads:
        2
    log:
        log_file=outdir+"/04_hisat2_mapped/log/"+"{sample}"+"_filt"+"_hisat2_mapped.log"
    params:
        shell_script=outdir+"/04_hisat2_mapped/log/"+"{sample}"+"_filt"+"_hisat2_mapped.sh",
        threads=config["hisat2_para"]["threads"],
        software=config["software"]["hisat2"],
        index=config["hisat2_para"]["index"],
        
    shell:
        """
        echo -e '{params.software} -x {params.index} -1 {input.filt_R1} -2 {input.filt_R2} -S {output.sam_file} -p {params.threads} -t' >& {params.shell_script}
        
        {params.software} -x {params.index} -1 {input.filt_R1} -2 {input.filt_R2} -S {output.sam_file} -p {params.threads} -t >& {log.log_file} 
        
        """