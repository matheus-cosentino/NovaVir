# workflow/rules/spades.smk

import os
from glob import glob

# --- Helper Functions ---

# This helper function still correctly defines the file dependencies for SPAdes.
def get_spades_input(wildcards):
    """
    Reads the layout file from the 'determine_sra_layout' checkpoint 
    and returns the appropriate dictionary of input files for SPAdes.
    """
    # This function is only called if spades is selected, but we still need 
    # the checkpoint logic to function. (Assuming checkpoint determine_sra_layout 
    # exists in get_data.smk)
    try:
        layout_file = checkpoints.determine_sra_layout.get(sample=wildcards.sample).output[0]
        
        with open(layout_file) as f:
            layout = f.read().strip()
    except:
        # If the checkpoint output is not available (e.g., in FASTQ mode where it might not run)
        # We assume paired-end default for local FASTQ data.
        # This part of the logic may need fine-tuning depending on your FASTQ/SRA rules.
        layout = "PE" 

    r1_path = f"{config['output_dir']}/{wildcards.sample}/trimmed/{wildcards.sample}_1.fastq.gz"
    r2_path = f"{config['output_dir']}/{wildcards.sample}/trimmed/{wildcards.sample}_2.fastq.gz"
    orphans_path = f"{config['output_dir']}/{wildcards.sample}/trimmed/{wildcards.sample}_orphans.fastq.gz"

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

def get_final_contigs_input(wildcards):
    """
    Determines the input for the select_final_contigs rule based on the sample's mode.
    """
    mode = SAMPLE_MODES[wildcards.sample]
    if mode == 'CONTIGS':
        return f"{config['output_dir']}/{wildcards.sample}/assembly/{wildcards.sample}/linked_contigs.fasta"
    elif mode in ['SRA', 'FASTQ']:
        return f"{config['output_dir']}/{wildcards.sample}/assembly/{wildcards.sample}/spades_contigs.fasta"
    else:
        raise ValueError(f"Unknown sample mode: {mode}")

def get_spades_rule_inputs(wildcards):
    """
    Returns a dictionary of named inputs for the spades rule.
    Returns an empty dictionary if the rule should be skipped (for CONTIGS mode),
    relying on the DAG to not require the output for these samples.
    """
    if SAMPLE_MODES[wildcards.sample] == 'CONTIGS':
        # This rule's output is not required for CONTIGS samples, so it won't be run.
        # Returning an empty dict is fine.
        return {}
    
    inputs = get_spades_input(wildcards) # returns {'r1': ..., 'r2': ...}
    inputs['layout'] = f"{config['output_dir']}/{wildcards.sample}/sra_layout/{wildcards.sample}.layout.txt"
    return inputs


rule spades:
    input:
        unpack(get_spades_rule_inputs)
    output:
        # Use a unique, prefixed output name
        contigs=f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/spades_contigs.fasta"
    
    shadow: 
        "minimal" 
    params:
        outdir=f"{config['output_dir']}/{{sample}}/assembly/",
        extra=config["params"]["spades"]["extra"]
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_spades.log"
    conda:
        SPADES
    shell:
        """
        # Read the layout from the input file
        # Note: If layout is not required (e.g., in some FASTQ paths), 
        # this line may need to be wrapped in an if block. 
        # Assuming layout is always produced for SRA/FASTQ for now.
        LAYOUT=$(cat {input.layout})

        # Build the reads command string based on the layout
        if [ "$LAYOUT" = "PE" ]; then
            READ_ARGS="--pe1-1 {input.r1} --pe1-2 {input.r2} --pe1-s {input.orphans}"
        else
            READ_ARGS="--s1 {input.r1}"
        fi
        
        spades.py \
            {params.extra} \
            $READ_ARGS \
            -t {resources.threads} \
            -m {resources.mem_mb} \
            --only-assembler \
            -o {params.outdir} \
            &> {log}
        """


rule link_preassembled_contigs:
    message:
        "MODE=CONTIGS detected. Linking pre-assembled contigs from CONTIGS_DIR for sample {wildcards.sample}."
    output:
        # Use a unique, prefixed output name
        contigs=f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/linked_contigs.fasta"
    
    params:
        mode = lambda wildcards: SAMPLE_MODES[wildcards.sample]

    run:
        if params.mode != 'CONTIGS':
            raise Exception(f"Rule link_preassembled_contigs called for sample {wildcards.sample} with mode {params.mode}, but it only supports CONTIGS mode.")
        
        # Caminho completo do arquivo de entrada FASTA/FNA
        try:
            # Note: CONTIGS_DIR must be globally available (defined in the main Snakefile)
            input_glob = os.path.join(CONTIGS_DIR, f"{wildcards.sample}*")
            possible_files = glob(input_glob + ".fasta") + glob(input_glob + ".fna")
            if not possible_files:
                raise IndexError
            input_file = possible_files[0]
        except IndexError:
            # We don't need to check SAMPLE_MODES here, as the input check handles skipping.
            # We only need to check for the file's existence if we are scheduled.
            raise FileNotFoundError(f"Error: Could not find pre-assembled contig file for {wildcards.sample} in {CONTIGS_DIR}/ (searched for: {wildcards.sample}*.fasta/fna)")

        # 1. Cria o diretório de saída esperado 
        assembly_dir = os.path.dirname(output.contigs)
        os.makedirs(assembly_dir, exist_ok=True)
        
        # 2. Cria o link simbólico
        shell(f"ln -sf $(readlink -f {input_file}) {output.contigs}")


rule select_final_contigs:
    message:
        "Selecting final contigs.fasta from mode-specific assembly file for {wildcards.sample}."
    
    output:
        f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/contigs.fasta"
    input:
        get_final_contigs_input
    
    run:
        source_file = input[0]
        # Create a symbolic link from the mode-specific output to the generic output
        print(f"Linking {os.path.basename(source_file)} to {os.path.basename(output[0])} for {wildcards.sample}...")
        shell(f"ln -sf $(readlink -f {source_file}) {output[0]}")