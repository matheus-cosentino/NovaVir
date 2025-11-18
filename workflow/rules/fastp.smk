rule fastp:
    input:
        r1="raw_data/{sample}_1.fastq.gz",
        r2="raw_data/{sample}_2.fastq.gz",
        layout=f"{config['output_dir']}/{{sample}}/sra_layout/{{sample}}.layout.txt"
    output:
        # Define all possible output files explicitly and statically.
        r1=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}_2.fastq.gz",
        orphans=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}_orphans.fastq.gz",
        html=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.html",
        json=f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.json"
    shadow: 
        "minimal" 
    log:
        f"{config['output_dir']}/{{sample}}/{{sample}}_fastp.log"
    params:
        # We only need to define non-file parameters here now.
        length_required=config['params']['fastp']['length_required'],
        quality=config['params']['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        LAYOUT=$(cat {input.layout})

        # Create empty placeholder files first to ensure all outputs exist.
        # This is the key to satisfying Snakemake for the single-end case.
        touch {output.r2} {output.orphans}

        # Build the command arguments based on the layout
        if [ "$LAYOUT" = "PE" ]; then
            CMD_ARGS="--in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --unpaired1 {output.orphans} --unpaired2 {output.orphans}"
        else
            CMD_ARGS="--in1 {input.r1} --out1 {output.r1}"
        fi

        # Run fastp. Note there is NO manual log redirection (&>).
        fastp \
            $CMD_ARGS \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality}
        """
