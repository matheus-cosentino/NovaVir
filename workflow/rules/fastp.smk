# workflow/rules/fastp.smk

#priority of rules
ruleorder: fastp_paired > fastp_unpaired

#rule to process paired data
rule fastp_paired:
    input:
        r1 = get_input_r1,
        r2 = get_input_r2
    output:
        # Define all possible output files explicitly and statically.
        r1= os.path.join(OUT_DIR, "{sample}/trimmed/{sample}_1.fastq.gz"),
        r2= os.path.join(OUT_DIR, "{sample}/trimmed/{sample}_2.fastq.gz"),
        orphans=os.path.join(OUT_DIR, "{sample}/trimmed/{sample}_orphans.fastq.gz")
        html=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.html",
        json=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.json"

    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_fastp.log"
    params:
        length_required=config['params']['fastp']['length_required'],
        quality=config['params']['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --unpaired1 {output.orphans} --unpaired2 {output.orphans} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality} \
            2> {log}
        """

#rule to process unpaired data
rule fastp_unpaired:
    input:
        reads = get_input_r1,
    output:
        r1=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}_1.fastq.gz",
        html=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.html",
        json=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.json"
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_fastp.log"
    params:
        length_required=config['params']['fastp']['length_required'],
        quality=config['params']['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        fastp \
            --in1 {input.reads} --out1 {output.reads} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality} \
            2> {log}
        """
