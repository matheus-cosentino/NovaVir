###################################################################################
#                       workflow/rules/get_data.smk                               #
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

###################
# --- HELPERS --- #
###################

PAIRED_SRA = [s for s, m in SAMPLE_META.items() if m['mode'] == 'SRA' and len(m['files']) == 2]
SINGLE_SRA = [s for s, m in SAMPLE_META.items() if m['mode'] == 'SRA' and len(m['files']) == 1]

##############################################
# --- 1. Download SRA Data (PAIRED) --- #
############################################# 

rule download_sra_data_paired:
    output:
        r1 = os.path.join(config["data_dir"], "{sample}_1.fastq.gz"),
        r2 = os.path.join(config["data_dir"], "{sample}_2.fastq.gz")
    log:
        os.path.join(config["output_dir"], "{sample}", "log", "{sample}_download.log")
    conda:
        DOWNLOAD
    shadow: 
        "minimal" 
    params:
        out_dir = config["data_dir"]
    wildcard_constraints:
        sample = "|".join(PAIRED_SRA) if PAIRED_SRA else "NO_PAIRED_SAMPLES"
    shell:
        """
        echo "Starting PAIRED download (Low Disk Mode) for {wildcards.sample}..." > {log}
        
        # MUDANÇA CRÍTICA:
        # Usamos 'fastq-dump' (antigo) em vez de 'fasterq-dump'.
        # --gzip: Escreve JÁ COMPRIMIDO (economiza muito espaço)
        # --split-3: Separa R1 e R2 corretamente
        
        fastq-dump \
            --split-3 \
            --gzip \
            --outdir {params.out_dir} \
            {wildcards.sample} >> {log} 2>&1
        """

##############################################
# --- 2. Download SRA Data (SINGLE) --- #
############################################# 

rule download_sra_single:
    output:
        r1 = os.path.join(config["data_dir"], "{sample}_1.fastq.gz")
    log:
        os.path.join(config["output_dir"], "{sample}", "log", "{sample}_download_single.log")
    conda:
        DOWNLOAD
    wildcard_constraints:
        sample = "|".join(SINGLE_SRA) if SINGLE_SRA else "NO_SINGLE_SAMPLES"
    shadow: 
        "minimal"
    params:
        out_dir = config["data_dir"]
    shell:
        """
        echo "Starting SINGLE download (Low Disk Mode) for {wildcards.sample}..." > {log}
        
        fastq-dump \
            --gzip \
            --outdir {params.out_dir} \
            {wildcards.sample} >> {log} 2>&1
            
        # Renomeia se necessário para garantir o padrão _1.fastq.gz se o fastq-dump gerar sem sufixo
        # O fastq-dump single as vezes gera apenas sample.fastq.gz
        
        if [ -f "{params.out_dir}/{wildcards.sample}.fastq.gz" ]; then
            mv "{params.out_dir}/{wildcards.sample}.fastq.gz" {output.r1}
        fi
        """