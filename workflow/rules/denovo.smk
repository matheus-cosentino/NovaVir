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
#                              version: 03.2026                                   #
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
        contigs = os.path.join(OUT_DIR, "{sample}", "spades", "kmer_{kmer_val}", "contigs.fasta")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "kmer_{kmer_val}_spades_assembly.log")
    params:
        outdir = os.path.join(OUT_DIR, "{sample}", "spades", "kmer_{kmer_val}"),
        input_args = get_spades_params,
        extra = config["spades"]["algorithm"],
        #kmer = lambda wildcards: wildcards.kmer_val
        kmer = lambda wildcards: wildcards.kmer_val.replace("_", ",")
    threads: 1
    conda:
        SPADES_ENV
    shadow:
        "minimal"
    shell:
        """
        # 1. Define Local Temp Directory
        # This uses the node's local disk (fast & no network lag)
        # If SLURM_TMPDIR isn't defined, it falls back to /tmp
        #define memory use within spades 
        mem_gb=$(({resources.mem_mb} / 1024))

        # 2. Run SPAdes with --tmp-dir pointing to LOCAL storage
        spades.py \
            {params.extra} \
            --threads {threads} \
            --memory $mem_gb \
            {params.input_args} \
            -o {params.outdir} \
            -k {params.kmer} \
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
        contigs = os.path.join(OUT_DIR, "{sample}", "megahit", "final.contigs.fa")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "meghit_assembly.log")
    params:
        outdir = os.path.join(OUT_DIR, "{sample}", "megahit"),
        input_args = get_megahit_params,
    threads: 1
    conda:
        DENOVO
    shell:
        """
        megahit \
          {params.input_args} \
          -o {params.outdir} \
          -t {threads} \
          -f \
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
        contigs = os.path.join(OUT_DIR, "{sample}", "flye", "assembly.fasta")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "flye_assembly.log")
    params:
        outdir = os.path.join(OUT_DIR, "{sample}", "flye"),
        extra = config["flye"]["type"]
    threads: 1
    conda:
        DENOVO
    shell:
        """
        # Flye fails if output directory exists
        rm -rf {params.outdir}

        # Run commands 

        flye --{params.extra} {input.reads} \
         --out-dir {params.outdir} \
         --threads {threads} \
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
        contigs = os.path.join(OUT_DIR, "{sample}", "raven", "assembly.fasta")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "raven_assembly.log")
    threads: 1
    conda:
        DENOVO
    shell:
        """
        # Run commands 
        raven --threads {threads} {input.reads} > {output.contigs} 2> {log}
        """

###############################
# --- 5. Medaka Correcton --- #
###############################

rule medaka_polish:
    input:
        reads = get_ONP_input,
        draft = os.path.join(OUT_DIR, "{sample}", "{assembler}", "assembly.fasta")
    output:
        consensus = os.path.join(OUT_DIR, "{sample}", "medaka_{assembler}", "consensus.fasta")
    params:
        # Defina o modelo aqui ou no config.yaml
        model = config["medaka_model"],
        outdir = os.path.join(OUT_DIR, "{sample}", "medaka_{assembler}")
    threads: 1
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "medaka_{assembler}.log")
    shell:
        """
        rm -rf {params.outdir}
        medaka_consensus -i {input.reads} \
                         -d {input.draft} \
                         -o {params.outdir} \
                         -t {threads} \
                         -m {params.model} \
        > {log} 2>&1
        
        """