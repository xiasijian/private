########################### fastp to filter data #################################
rule fastp_filt:
    input:
        raw_R1=data_dir+"{sample}"+"/"+"{sample}"+"_1"+fastq_suffix,
        raw_R2=data_dir+"{sample}"+"/"+"{sample}"+"_2"+fastq_suffix
    output:
        filt_R1=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R1_filt"+".fq.gz",
        filt_R2=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R2_filt"+".fq.gz",
        json=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_filt"+"_fastp.json"
    threads:8

    log:
        log_file=outdir+"/02_fastp_filt/"+"log/"+"{sample}"+"fastp_filt"+".log"
    params:
        shell_script=outdir+"/02_fastp_filt/"+"log/"+"{sample}"+"_fastp_filt"+".sh",
        quality=config["fastp_para"]["quality"],
        software=config["software"]["fastp"],
        filt_read_length=config["fastp_para"]["filt_read_length"],
        thread=config["fastp_para"]["thread"],
        html = outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_filt"+"_fastp.html",
        fastp_env=config["software"]["fastp_env"]
    shell:
        """
        ## source activate {params.fastp_env}
        ## {params.software} --length_required=45 --trim_poly_g --cut_tail --correction -w {params.thread} -h {params.html} -j {output.json} -i {input.raw_R1} -o {output.filt_R1} -I {input.raw_R2} -O {output.filt_R2} 2> {log.log_file}
        echo -e "source activate {params.fastp_env} \n{params.software} --length_required=45 --trim_poly_g --cut_tail --correction -w {params.thread} -h {params.html} -j {output.json} -i {input.raw_R1} -o {output.filt_R1} -I {input.raw_R2} -O {output.filt_R2}" > {params.shell_script}
       	sh {params.shell_script} > {log.log_file}
	"""
