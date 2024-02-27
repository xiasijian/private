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
        r1_file = DATA_DIR + "{experiment}"+"/" +  "{experiment}" + "_1.fastq.gz",
        r2_file = DATA_DIR + "{experiment}"+"/" +  "{experiment}" + "_2.fastq.gz",
        start_file = outdir + "/"+"{experiment}" + ".txt"
    output:
        r1_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_1.fq.gz",
        r2_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_2.fq.gz",
        html = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_fastp.html",
        json = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_fastp.json",
        clean_fastq_completed = outdir + "/01_clean_fastq/output/" + "{experiment}" +"_completed.txt"
    log:
        log_file = outdir + "/01_clean_fastq/"+"log/"+ "{experiment}" + "_fastp_clean.log"
    params:
        mcc_env =config["conda_env"]["mcc_env"],
	    fastp =config["software"]["fastp"],
        shell_script =outdir + "/01_clean_fastq/"+"log/"+"{experiment}" + "_fastp_clean.sh",
        threads =config["fastp_para"]["threads"],
        quality_cutoff = config["fastp_para"]["quality"],
        minReadLength = config["fastp_para"]["filt_read_length"]
    threads:
        config["fastp_para"]["threads"]
    shell:
        """
        echo -e '{params.fastp} -q {params.quality_cutoff} -w {params.threads} -i {input.r1_file} -I {input.r2_file} -o {output.r1_filt} -O {output.r2_filt} -h {output.html} -j {output.json}' > {params.shell_script}

        bash {params.shell_script} >& {log.log_file}

        touch {output.clean_fastq_completed}

        """


# rule clean_fastq:
#     input:
#         # r1_file = lambda wildcards: DATA_DIR + get_fq1_from_experiment[wildcards.experiment],
#         # r2_file = lambda wildcards: DATA_DIR + get_fq2_from_experiment[wildcards.experiment],
#         r1_file = DATA_DIR + "{experiment}"+"/" +  "{experiment}" + "_1.fastq.gz",
#         r2_file = DATA_DIR + "{experiment}"+"/" +  "{experiment}" + "_2.fastq.gz",
#         start_file = outdir + "/"+"{experiment}" + ".txt"
#     output:
#         r1_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_1.fq.gz",
#         r2_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_2.fq.gz",
#         clean_fastq_completed = outdir + "/01_clean_fastq/output/" + "{experiment}" +"_completed.txt"
#     log:
#         log_file = outdir + "/01_clean_fastq/"+"log/"+ "{experiment}" + "_trim_galore_clean.log"
#     params:
#         mcc_env =config["conda_env"]["mcc_env"],
# 	    trim_galore =config["software"]["trim_galore"],
#         shell_script =outdir + "/01_clean_fastq/"+"log/"+"{experiment}" + "_trim_galore_clean.sh",
#         threads =config["trim_galore_para"]["threads"],
#         quality_cutoff = config["trim_galore_para"]["quality_cutoff"],
#         minReadLength = config["trim_galore_para"]["minReadLength"],
# 	    res_dir=outdir+"/01_clean_fastq/"+"output/",
#         basename="{experiment}",
#         barcode_para = barcode_command
#     threads:
#         config["trim_galore_para"]["threads"]
#     shell:
#         """
#         echo -e 'source activate {params.mcc_env} \n{params.trim_galore} -j {params.threads} --output_dir {params.res_dir}{params.barcode_para} -q {params.quality_cutoff} --length {params.minReadLength} --paired {input.r1_file} {input.r2_file} --basename {params.basename} --no_report_file' > {params.shell_script}

#         bash {params.shell_script} >& {log.log_file}

#         touch {output.clean_fastq_completed}

#         """
    
# rule fastqc_by_fastp:
#     input:
#         r1_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_1.fq.gz",
#         r2_filt = outdir + "/01_clean_fastq/output/" + "{experiment}" + "_val_2.fq.gz",
#         clean_fastq_completed = outdir + "/01_clean_fastq/output/" + "{experiment}" +"_completed.txt"
#     output:
#         fastp_html = outdir+"/01_clean_fastq/"+"output/"+"{experiment}.filt.fastp.qc.html"
#     log:
#         outdir+"/01_clean_fastq/"+"log/" + "{experiment}"+".filt.fastp.qc.log"
#     params:
#         mcc_env =config["conda_env"]["mcc_env"],
#         chipseq_env = config["conda_env"]["chipseq_env"],
#         fastp = config["software"]["fastp"],
#         shell_script = outdir+"/01_clean_fastq/"+"log/" + "{experiment}"+".filt.fastp.qc.sh"
#     shell:
#         """
#         echo -e 'source activate {params.chipseq_env} \n{params.fastp} -i {input.r1_filt} -I {input.r2_filt} -h {output.fastp_html}' > {params.shell_script}

#         bash {params.shell_script} >& {log}

#         """
        
