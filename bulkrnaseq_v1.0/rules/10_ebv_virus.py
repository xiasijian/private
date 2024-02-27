rule venus:
    input:
        filt_R1=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R1_filt"+".fq.gz",
        filt_R2=outdir+"/02_fastp_filt/"+"output/"+"{sample}"+"_R2_filt"+".fq.gz"
    output:
        venus_out=outdir+"/10_venus/"+"{sample}"+"/detection/bulk_paired-end/detection_output.tsv"
    threads: 5
    params:
        venus_env=config["conda_env"]["venus_env"],
        venus_script="/share/home/xiasj/01_biosoft/Venus/src/module-detection/module-detection.py",
        shell_script=outdir+"/10_venus/"+"{sample}"+"_run_venus.sh",
        virus_chr="/share/home/xiasj/01_biosoft/Venus/reference_files/ebv_chr_ref.tsv",
        virusGenome="/share/home/xiasj/01_biosoft/Venus/indices/ebv_NC_007605.genomeDir",
        humanGenome="/share/home/xiasj/01_biosoft/Venus/indices/human2.genomeDir",
        outdir=outdir+"/10_venus/"+"{sample}"+"/detection/bulk_paired-end/"
    shell:
        """
        echo -e 'source activate {params.venus_env} \npython3 {params.venus_script} --read {input.filt_R1} {input.filt_R2} --virusThreshold 5 --virusChrRef {params.virus_chr} --virusGenome {params.virusGenome} --humanGenome {params.humanGenome} --readFilesCommand zcat --thread 5 --out {params.outdir}' > {params.shell_script}
        
        sh {params.shell_script}
        
        """