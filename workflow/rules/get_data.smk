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
#                              version: 03.2026                                   #
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
        out_dir = config["data_dir"],
        # FIX 1: Corrected Syntax with f-string
        tmpdir = lambda wildcards: os.path.join(config["data_dir"], f"tmp_{wildcards.sample}")
    wildcard_constraints:
        sample = "|".join(PAIRED_SRA) if PAIRED_SRA else "NO_PAIRED_SAMPLES"
    shell:
        """
        exec &> {log}
        echo "Step 1: Prefetching SRA file for {wildcards.sample}..."
        rm -rf ~/ncbi/public/sra/{wildcards.sample}.sra*
        prefetch {wildcards.sample} --max-size 100G
        
        echo "Step 2: Validating SRA integrity..."
        vdb-validate {wildcards.sample}

        echo "Step 3: Running fasterq-dump from local SRA file..."
        mkdir -p {params.tmpdir}
        fasterq-dump {wildcards.sample} \
            --split-files \
            --threads {threads} \
            --temp {params.tmpdir} \
            --outdir {params.tmpdir} \
            --skip-technical

        echo "Step 4: Compressing with pigz (faster) or gzip..."
        # Se 'pigz' estiver disponível no cluster, use-o para ganhar tempo
        pigz -p {threads} "{params.tmpdir}/{wildcards.sample}_1.fastq" || gzip "{params.tmpdir}/{wildcards.sample}_1.fastq"
        pigz -p {threads} "{params.tmpdir}/{wildcards.sample}_2.fastq" || gzip "{params.tmpdir}/{wildcards.sample}_2.fastq"
        
        echo "Step 5: Moving files and cleaning up..."
        mv "{params.tmpdir}/{wildcards.sample}_1.fastq.gz" {output.r1}
        mv "{params.tmpdir}/{wildcards.sample}_2.fastq.gz" {output.r2}
        
        # Remove o arquivo .sra baixado para economizar espaço no diretório do NCBI
        rm -rf ~/ncbi/public/sra/{wildcards.sample}.sra
        rm -rf {params.tmpdir}
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
    # FIX 4: Added missing params section
    params:
        out_dir = config["data_dir"],
        tmpdir = lambda wildcards: os.path.join(config["data_dir"], f"tmp_{wildcards.sample}")
    shell:
        """
        echo "Creating temp dir for download {wildcards.sample}..." > {log}
        mkdir -p {params.tmpdir}

        echo "Starting SingleEnd download for {wildcards.sample}..." >> {log}
        fasterq-dump --split-files --threads {resources.threads} -O {params.tmpdir} {wildcards.sample} >> {log} 2>&1

        echo "Compressing and renaming..." >> {log}
        
        # Logic to handle naming variations (sample.fastq vs sample_1.fastq)
        # Note: We must compress BEFORE moving to .gz output
        
        if [ -f "{params.tmpdir}/{wildcards.sample}_1.fastq" ]; then
            gzip "{params.tmpdir}/{wildcards.sample}_1.fastq"
            mv "{params.tmpdir}/{wildcards.sample}_1.fastq.gz" {output.r1}
        elif [ -f "{params.tmpdir}/{wildcards.sample}.fastq" ]; then
            gzip "{params.tmpdir}/{wildcards.sample}.fastq"
            mv "{params.tmpdir}/{wildcards.sample}.fastq.gz" {output.r1}
        fi

        echo "Cleaning up temp dir..." >> {log}
        rm -rf {params.tmpdir}
        """