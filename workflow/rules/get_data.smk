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

PAIRED_SRA_LIST = [s for s, m in SAMPLE_META.items() if m['mode'] == 'SRA' and len(m['files']) == 2]
SINGLE_SRA_LIST = [s for s, m in SAMPLE_META.items() if m['mode'] == 'SRA' and len(m['files']) == 1]

##############################################
# --- 1. Download SRA Data of Libraries --- #
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
        #r2 = os.path.join(config["data_dir"], "{sample}_2.fastq.gz")
    wildcard_constraints:
        sample = "|".join(PAIRED_SRA) if PAIRED_SRA else "NO_PAIRED_SAMPLES"
    shadow: 
        "minimal"
    shell:
         """
        echo "Starting PAIRED download for {wildcards.sample}..." > {log}
        
        fasterq-dump --split-files --threads {threads} -O . {wildcards.sample} >> {log} 2>&1

        # Compress
        gzip "{wildcards.sample}_1.fastq"
        gzip "{wildcards.sample}_2.fastq"
        
        # Move to final output
        mv "{wildcards.sample}_1.fastq.gz" {output.r1}
        mv "{wildcards.sample}_2.fastq.gz" {output.r2}
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
    # This constraint ensures this rule ONLY runs for samples we identified as Single
    wildcard_constraints:
        sample = "|".join(SINGLE_SRA) if SINGLE_SRA else "NO_SINGLE_SAMPLES"
    shadow: 
        "minimal"
    shell:
        """
        echo "Starting SINGLE download for {wildcards.sample}..." > {log}
        
        fasterq-dump --split-files --threads {threads} -O . {wildcards.sample} >> {log} 2>&1

        # fasterq-dump outputs just 'sample.fastq' for single end, OR 'sample_1.fastq' depending on version.
        # We handle both cases safely:
        
        if [ -f "{wildcards.sample}_1.fastq" ]; then
            gzip "{wildcards.sample}_1.fastq"
            mv "{wildcards.sample}_1.fastq.gz" {output.r1}
        elif [ -f "{wildcards.sample}.fastq" ]; then
            gzip "{wildcards.sample}.fastq"
            mv "{wildcards.sample}.fastq.gz" {output.r1}
        fi
        """