# workflow/rules/spades.smk

# This helper function still correctly defines the file dependencies for Snakemake's DAG.
def get_spades_input(wildcards):
    """
    Reads the layout file from the 'determine_sra_layout' checkpoint 
    and returns the appropriate dictionary of input files for SPAdes.
    """
    layout_file = checkpoints.determine_sra_layout.get(sample=wildcards.sample).output[0]
    
    with open(layout_file) as f:
        layout = f.read().strip()

    r1_path = f"{config['output_dir']}/{{sample}}/trimmed/{wildcards.sample}_1.fastq.gz"
    r2_path = f"{config['output_dir']}/{{sample}}/trimmed/{wildcards.sample}_2.fastq.gz"
    orphans_path = f"{config['output_dir']}/{{sample}}/trimmed/{wildcards.sample}_orphans.fastq.gz"

    if layout == "PE":
        return {
            "r1": r1_path,
            "r2": r2_path,
            "orphans": orphans_path
        }
    else: # SE
        return {
            "r1": r1_path
        }

# This single rule now correctly handles both PE and SE assembly.
rule spades:
    input:
        # Unpack the required read files based on the layout
        unpack(get_spades_input),
        # **Explicitly add the layout file as an input to use in the shell command**
        layout=f"{config['output_dir']}/{{sample}}/sra_layout/{{sample}}.layout.txt"
    output:
        contigs=f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/contigs.fasta"
    shadow: 
        "minimal" 
    params:
        outdir=f"{config['output_dir']}/{{sample}}/assembly/",
        # The 'extra' parameter is fine as it is
        extra=config["params"]["spades"]["extra"]
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_spades.log"
    conda:
        SPADES
    shell:
        """
        # Read the layout from the input file
        LAYOUT=$(cat {input.layout})

        # Build the reads command string based on the layout
        if [ "$LAYOUT" = "PE" ]; then
            READ_ARGS="--pe1-1 {input.r1} --pe1-2 {input.r2} --pe1-s {input.orphans}"
        else
            READ_ARGS="--s1 {input.r1}"
        fi

        # Note: SPAdes has updated its read input flags.
        # -1, -2, -s are for single libraries. Using --pe<#>-<1,2,s> is more robust.
        spades.py \
            {params.extra} \
            $READ_ARGS \
            -t {resources.threads} \
            -m {resources.mem_mb} \
            --only-assembler \
            -o {params.outdir} \
            &> {log}
        """

