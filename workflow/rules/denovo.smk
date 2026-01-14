###################################################################################
#                       workflow/rules/denovo.smk                                 #
#                         MSc. Matheus Cosentino                                  # 
###################################################################################
#                                                                                 #
# oooooooooo.    o8o                               oooooo     oooo  o8o           #
# `888'   `Y8b   `"'                                `888.     .8'   `"'           #
#  888      888 oooo   .oooo.o  .ooooo.   .ooooo.    `888.   .8'   oooo  oooo d8b #
#  888      888 `888  d88(  "8 d88' `"Y8 d88' `88b    `888. .8'    `888  `888""8P #
#  888      888  888  `"Y88b.  888       888   888     `888.8'      888   888     #
#  888     d88'  888  o.  )88b 888   .o8 888   888      `888'       888   888     #
# o888bood8P'   o888o 8""888P' `Y8bod8P' `Y8bod8P'       `8'       o888o d888b    #
#                                                                                 #
###################################################################################
#                              version: 12.2025                                   #
###################################################################################

###############################
# --- 1. Spades Assembly --- #
##############################

rule spades:
    input:
        r1 = get_denovo_r1,
        r2 = get_denovo_r2,
        extra = get_denovo_unpaired
    output:
        contigs = "{out_dir}/{sample}/spades/kmer_{kmer_val}/contigs.fasta"
    log:
        "{out_dir}/{sample}/log/kmer_{kmer_val}_spades_assembly.log"
    params:
        outdir = "{out_dir}/{sample}/spades/kmer_{kmer_val}",
        input_args = get_spades_params,
        extra = config["spades"]["algorithm"],
        kmer = lambda wildcards: wildcards.kmer_val
    conda:
        DENOVO
    shadow:
        "minimal"
    shell:
        """
        # 1. Define Local Temp Directory
        # This uses the node's local disk (fast & no network lag)
        # If SLURM_TMPDIR isn't defined, it falls back to /tmp
        LOCAL_TMP=${{SLURM_TMPDIR:-/tmp}}
        
        #where spades temp data will be written
        echo "[INFO] SPAdes Local TMP is set to: $LOCAL_TMP" > {log}
        
        # Check how much space is available there (crucial for SPAdes)
        df -h "$LOCAL_TMP" >> {log}
        
        #define memory use within spades 
        mem_gb=$(({resources.mem_mb} / 1024))

        # 2. Run SPAdes with --tmp-dir pointing to LOCAL storage
        spades.py \
            {params.extra} \
            --threads {resources.threads} \
            --memory $mem_gb \
            {params.input_args} \
            -o {params.outdir} \
            -k {params.kmer} \
            --tmp-dir "$LOCAL_TMP" \
            --only-assembler \
            > {log} 2>&1
        """

###############################
# --- 2. Megahit Assembly --- #
##############################

rule megahit:
    input:
        r1 = get_denovo_r1,
        r2 = get_denovo_r2,
        extra = get_denovo_unpaired
    output:
        contigs = "{out_dir}/{sample}/megahit/final.contigs.fa"
    log:
        "{out_dir}/{sample}/log/meghit_assembly.log"
    params:
        outdir = "{out_dir}/{sample}/megahit",
        input_args = get_megahit_params,
    conda:
        DENOVO
    shell:
        """
        megahit \
          {params.input_args} \
          -o {params.outdir} \
          -t {resources.threads} \
          > {log} 2>&1
             
        """


#############################
# --- 3. Flye Assembly --- #
############################

rule flye:  
    """
    WARNING: Flye is designed for Long Reads (Nanopore/PacBio).
    Do not use Illumina reads here.
    """
    input:
        reads = get_ONP_input
    output:
        contigs = "{out_dir}/{sample}/flye/assembly.fasta"
    log:
        "{out_dir}/{sample}/log/flye_assembly.log"
    params:
        outdir = "{out_dir}/{sample}/flye",
        extra = config["flye"]["type"]
    conda:
        DENOVO
    shell:
        """
        # Flye fails if output directory exists
        rm -rf {params.outdir}

        # Run commands 

        flye --{params.extra} {input.reads} \
         --out-dir {params.outdir} \
         --threads {resources.threads} \
         --meta \
         > {log} 2>&1
        """

#############################
# --- 4. Raven Assembly --- #
############################

rule raven:  
    """
    WARNING: Flye is designed for Long Reads (Nanopore/PacBio).
    Do not use Illumina reads here.
    """
    input:
        reads = get_ONP_input
    output:
        contigs = "{out_dir}/{sample}/raven/assembly.fasta"
    log:
        "{out_dir}/{sample}/log/raven_assembly.log"
    conda:
        DENOVO
    shell:
        """
        # Run commands 
        raven --threads {resources.threads} {input.reads} > {output.contigs} 2> {log}
        """

###############################
# --- 5. Medaka Correcton --- #
###############################

rule medaka_polish:
    input:
        reads = get_ONP_input,
        draft = "{out_dir}/{sample}/{assembler}/assembly.fasta"
    output:
        consensus = "{out_dir}/{sample}/medaka_{assembler}/consensus.fasta"
    params:
        # Defina o modelo aqui ou no config.yaml
        model = config["medaka_model"],
        outdir = "{out_dir}/{sample}/medaka_{assembler}"
    log:
        "{out_dir}/{sample}/log/medaka_{assembler}.log"
    shell:
        """
        rm -rf {params.outdir}
        medaka_consensus -i {input.reads} \
                         -d {input.draft} \
                         -o {params.outdir} \
                         -t {resources.threads} \
                         -m {params.model} \
        > {log} 2>&1
        
        """