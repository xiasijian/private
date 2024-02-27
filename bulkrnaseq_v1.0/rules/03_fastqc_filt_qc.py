########################### fastqc for filter data #################################
rule fastqc_filt:
    input:
        filt_R1=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R1_filt"+".fq.gz",
        filt_R2=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R2_filt"+".fq.gz"
    output:
        fastqc_filt=outdir+"/03_fastqc_filt/"+"output/"+"{sample}"+"fastq_filt"+"_completed.txt"
    log:
        log_file=outdir+"/03_fastqc_filt/"+"log/"+"{sample}"+"fastq_filt"+"_completed.log"
    params:
        shell_script=outdir+"/03_fastqc_filt/"+"log/"+"{sample}"+"fastq_filt"+"_completed.sh",
        out_dir=outdir+"/03_fastqc_filt/"+"output/",
        software=config["software"]["fastqc"],
        threads=config["fastqc_para"]["threads"],
        outfile=outdir+"/03_fastqc_filt/output/"+"{sample}"+"fastq_filt"+"_completed.txt",
        java_env=config["software"]["java_env"]
    
    threads: 2
    
    shell:
        """
        echo -e 'source activate {params.java_env} \n{params.software} -f {input.filt_R1} {input.filt_R2} -o {params.out_dir}\ncd {params.out_dir}\ntouch {params.outfile}'>{params.shell_script}
        
        source activate {params.java_env}
        {params.software} -t {params.threads} -o {params.out_dir} {input.filt_R1} {input.filt_R2}
                
        touch {params.outfile}
        
        """