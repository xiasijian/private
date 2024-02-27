rule call_narrow_peaks:
    input:
        treatment =  lambda wildcards outdir + "/02_align/" + get_treat_from_sample[wildcards.sample] + "_bowtie2_mapped_sorted_rmdp.bam"
        control = lambda wildcards outdir + "/02_align/" + get_control_from_sample[wildcards.sample]  + "_bowtie2_mapped_sorted_rmdp.bam",
    log:
        outdir + "/03_macs2/" + "{sample}" + "_macs2_call_narrow_peaks.log"
    output:
        macs2_narrow_peaks = outdir + "/03_macs2/" + "{sample}" + "_peaks.narrowPeak"
    threads: 5
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        macs2 = config["software"]["macs2"],
        name = "{sample}"
        format = str(config['macs2_para']['format']),
        genomesize = str(config['macs2_para']['genomesize']),
        qvalue = str(config['macs2_para']['qvalue'])
        shell_script = outdir + "/03_macs2/" + "{sample}" + "_macs2_call_narrow_peaks.sh"
        
    shell:
    """
    echo -e '
    {params.macs2} callpeak -t {input.treatment} -c {input.control} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir {params.outdir}
    ' > {params.shell_script}
    
    sh {params.shell_script} >& {log}
    """
    
rule call_broad_peaks:
    input:
        treatment =  lambda wildcards outdir + "/02_align/" + get_treat_from_sample[wildcards.sample] + "_bowtie2_mapped_sorted_rmdp.bam",
        control = lambda wildcards outdir + "/02_align/" + get_control_from_sample[wildcards.sample]  + "_bowtie2_mapped_sorted_rmdp.bam"
    log:
        outdir + "/03_macs2/" + "{sample}" + "_macs2_call_broad_peaks.log"
    output:
        macs2_broad_peaks = outdir + "/03_macs2/" + "{sample}" + "_peaks.broadPeak"
    threads: 5
    params:
        chipseq_env = config["conda_env"]["chipseq_env"],
        macs2 = config["software"]["macs2"],
        name = "{sample}"
        format = str(config['macs2_para']['format']),
        genomesize = str(config['macs2_para']['genomesize']),
        qvalue = str(config['macs2_para']['qvalue'])
        shell_script = outdir + "/03_macs2/" + "{sample}" + "_macs2_call_broad_peaks.sh"
        
    shell:
    """
    echo -e '
    {params.macs2} callpeak -t {input.treatment} -c {input.control} {params.format} --broad --broad-cutoff 0.1 {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir {params.outdir}
    ' > {params.shell_script}
    
    sh {params.shell_script} >& {log}
    """