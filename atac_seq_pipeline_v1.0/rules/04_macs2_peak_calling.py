rule atacseq_peak_set1:
    input:
        #rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam",
        #rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam.bai"
        blacklist_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam",
        blacklist_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam.bai"
    output:
        macs2_peaks_set1 = outdir + "/03_macs2/" + "{experiment}" + "_set1_peaks.narrowPeak"
    log:
        outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        macs2 = config["software"]["macs2"],
        name = "{experiment}" + "_set1",
        outdir =  outdir + "/03_macs2/",
        format = str(config['macs2_para']['format']),
        genomesize = str(config['macs2_para']['genomesize']),
        qvalue = str(config['macs2_para']['qvalue']),
        set1 = str(config['macs2_para']['set1']["para"]),
        shell_script = outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks_set1.sh"
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'source activate {params.chipseq_env}' > {params.shell_script}
        echo -e '{params.macs2} callpeak -t {input.blacklist_bam} {params.format} \\
            {params.genomesize} --name {params.name} -q {params.qvalue} {params.set1} --outdir {params.outdir}' >> {params.shell_script}
        
        bash {params.shell_script}
        
        """

rule atacseq_peak_set2:
    input:
        #rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam",
        #rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam.bai"
        blacklist_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam",
        blacklist_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam.bai"
    output:
        macs2_peaks_set2 = outdir + "/03_macs2/" + "{experiment}" + "_set2_peaks.narrowPeak"
    log:
        outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        macs2 = config["software"]["macs2"],
        outdir =  outdir + "/03_macs2/",
        format = str(config['macs2_para']['format']),
        genomesize = str(config['macs2_para']['genomesize']),
        qvalue = str(config['macs2_para']['qvalue']),
        name = "{experiment}" + "_set2",
        set2 = str(config['macs2_para']['set2']["para"]),
        shell_script = outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks_set2.sh"
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'source activate {params.chipseq_env}' > {params.shell_script}
        echo -e '{params.macs2} callpeak -t {input.blacklist_bam} {params.format} \\
            {params.genomesize} --name {params.name} -q {params.qvalue} {params.set2} --outdir {params.outdir}' >> {params.shell_script}
        
        bash {params.shell_script}
        
        """

rule atacseq_peak_set3:
    input:
        #rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam",
        #rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam.bai"
        blacklist_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam",
        blacklist_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam.bai"
    output:
        macs2_peaks_set3 = outdir + "/03_macs2/" + "{experiment}" + "_set3_peaks.narrowPeak"
    log:
        outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks_set3.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        macs2 = config["software"]["macs2"],
        outdir =  outdir + "/03_macs2/",
        format = str(config['macs2_para']['format']),
        genomesize = str(config['macs2_para']['genomesize']),
        name = "{experiment}" + "_set3",
        set3 = str(config['macs2_para']['set3']["para"]),
        shell_script = outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks_set3.sh"
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'source activate {params.chipseq_env}' > {params.shell_script}
        echo -e '{params.macs2} callpeak -t {input.blacklist_bam} {params.format} \\
            {params.genomesize} --name {params.name} {params.set3} --outdir {params.outdir}' >> {params.shell_script}
        
        bash {params.shell_script} > {log}
        
        """

rule atacseq_peak_set5:
    input:
        #rmdp_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam",
        #rmdp_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp.bam.bai"
        blacklist_bam = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam",
        blacklist_bam_index = outdir + "/02_align/" + "{experiment}" + "_bowtie2_mapped_q30_rmdp_rmblacklist.bam.bai"
    output:
        macs2_peaks_set5 = outdir + "/03_macs2/" + "{experiment}" + "_set5_peaks.narrowPeak"
    log:
        outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks_set5.log"
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        macs2 = config["software"]["macs2"],
        outdir =  outdir + "/03_macs2/",
        name = "{experiment}" + "_set5",
        set5 = str(config['macs2_para']['set5']["para"]),
        shell_script = outdir + "/03_macs2/" + "{experiment}" + "_macs2_call_peaks_set5.sh"
    shell:
        """
        rm -f {params.shell_script}
        echo -e 'source activate {params.chipseq_env}' > {params.shell_script}
        echo -e '{params.macs2} callpeak -t {input.blacklist_bam} \\
            --name {params.name} {params.set5} --outdir {params.outdir}' >> {params.shell_script}
        
        bash {params.shell_script} > {log}
        
        """