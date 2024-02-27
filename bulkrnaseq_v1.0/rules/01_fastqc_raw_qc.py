rule fastqc_raw:
    input:
        raw_R1=data_dir+"{sample}"+"/"+"{sample}"+"_1"+fastq_suffix,
        raw_R2=data_dir+"{sample}"+"/"+"{sample}"+"_2"+fastq_suffix
    output:
        fastqc_raw=outdir+"/01_fastqc_raw/"+"output/"+"{sample}"+"fastq_raw_completed.txt"
    log:
        log_file=outdir+"/01_fastqc_raw/"+"log/"+"{sample}"+"fastq_raw_completed.log"
    params:
        shell_script=outdir+"/01_fastqc_raw/"+"log/"+"{sample}"+"fastq_raw_completed.sh",
        out_dir=outdir+"/01_fastqc_raw/"+"output/",
        software=config["software"]["fastqc"],
        outfile=outdir+"/01_fastqc_raw/"+"output/"+"{sample}"+"fastq_raw_completed.txt",
        threads=config["fastqc_para"]["threads"],
        java_env=config["software"]["java_env"]
    
    shell:
        """
        echo -e 'source activate {params.java_env} \n{params.software} -f {input.raw_R1} {input.raw_R2} -o {params.out_dir}\ncd {params.out_dir}\ntouch {params.outfile}'>{params.shell_script}
        
        source activate {params.java_env}
        {params.software} -t {params.threads} -o {params.out_dir} {input.raw_R1} {input.raw_R2} >& {log.log_file}
        
        touch {params.outfile}
        
        """