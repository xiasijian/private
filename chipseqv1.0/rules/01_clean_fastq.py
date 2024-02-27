adapter1 = config["trim_galore_para"]["adapter1"]
adapter2 = config["trim_galore_para"]["adapter2"]
if len(adapter1)!=0 and len(adapter2)==0:
    barcode_command = " -a " + adapter1
elif len(adapter1)==0 and len(adapter2)!=0:
    barcode_command = " -a " + adapter2
elif len(adapter1)==0 and len(adapter2)==0:
    barcode_command=""
else:
    barcode_command = " -a " + adapter1 + " " + "-a2 " + adapter2


print(barcode_command)


rule clean_fastq:
    input:
        r1_file = lambda wildcards: DATA_DIR + get_fq1_from_experiment[wildcards.experiment],
        r2_file = lambda wildcards: DATA_DIR + get_fq2_from_experiment[wildcards.experiment],
        start_file = outdir + "/"+"{experiment}" + ".txt"
    output:
        r1_filt = outdir + "/01_clean_fastq/" + "{experiment}" + "_val_1.fq.gz",
        r2_filt = outdir + "/01_clean_fastq/" + "{experiment}" + "_val_2.fq.gz"
    log:
        log_file = outdir + "/01_clean_fastq/"+"log/"+ "{experiment}" + "_trim_galore_clean.log"
    params:
        mcc_env =config["conda_env"]["mcc_env"],
	    trim_galore =config["software"]["trim_galore"],
        shell_script =outdir + "/01_clean_fastq/"+"log/"+"{experiment}" + "_trim_galore_clean.sh",
        threads =config["trim_galore_para"]["threads"],
        quality_cutoff = config["trim_galore_para"]["quality_cutoff"],
        minReadLength = config["trim_galore_para"]["minReadLength"],
	    res_dir=outdir+"/01_clean_fastq/"+"output/",
        basename="{experiment}",
        barcode_para = barcode_command
    threads:
        config["trim_galore_para"]["threads"]
    shell:
        """
        echo -e 'module load anaconda \nsource activate {params.mcc_env} \n{params.trim_galore} -j {params.threads} --output_dir {params.res_dir}{params.barcode_para} -q {params.quality_cutoff} --length {params.minReadLength} --paired {input.r1_file} {input.r2_file} --basename {params.basename}' > {params.shell_script}
        sh {params.shell_script} >& {log.log_file}

        """
    
rule fastqc_by_fastp:
    input:
        r1_filt = outdir + "/01_clean_fastq/" + "{experiment}" + "_val_1.fq.gz",
        r2_filt = outdir + "/01_clean_fastq/" + "{experiment}" + "_val_2.fq.gz"
    output:
        fastp_html = outdir+"/01_clean_fastq/"+"output/"+"{experiment}.filt.fastp.qc.html"
    log:
        outdir+"/01_clean_fastq/"+"output/" + "{experiment}.filt.fastp.qc.log"
    params:
        fastp = config["software"]["fastp"],
        shell_script = outdir+"/01_clean_fastq/"+"output/" + "{experiment}.filt.fastp.qc.sh"
    message:
        "Quality check of trimmed {wildcards.experiment} sample with fastp"
    shell:
    """
    echo -e 'module load anaconda \nsource activate {params.mcc_env} \n{params.fastp} -i {input.r1_filt} -I {input.r2_filt} -h {output.fastp_html}' > {params.shell_script}
    sh {params.shell_script}
    
    """
        