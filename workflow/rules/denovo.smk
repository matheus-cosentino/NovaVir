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
rule spades:
    input:
        r1 = get_denovo_r1,
        r2 = get_denovo_r2,
        extra = get_denovo_unpaired
    output:
        contigs = "{out_dir}/{sample}/spades/contigs.fasta"
    log:
        "{out_dir}/{sample}/log/spades_assembly.log"
    params:
        outdir = "{out_dir}/{sample}/spades",
        input_args = get_spades_params,
        extra = config["spades"]["algorithm"]
    conda:
        DENOVO
    shell:
        """
        
        mem_gb=$(({resources.mem_mb} / 1024))

        spades.py \
            {params.extra} \
            --threads {threads} \
            --memory $mem_gb \
            {params.input_args} \
            -o {params.outdir} \
            > {log} 2>&1
        """

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
        rm -rf {params.outdir}

        megahit.py \
          {params.input_args} \
          -o {params.outdir} \
          > {log} 2>&1
      
        """

rule flye:  
    """
    WARNING: Flye is designed for Long Reads (Nanopore/PacBio).
    Do not use Illumina reads here.
    """
    input:
        reads = get_flye_input
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

        flye \
         --{params.extra} {input.reads} \
         --out-dir {params.outdir} \
         --threads {resources.threads} \
         --meta \
         > {log} 2>&1
        """