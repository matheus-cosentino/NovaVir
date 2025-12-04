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


##############################################
# --- 1. Download SRA Data of Libraries --- #
############################################# 

rule download_sra_data:
    output:
        r1 = os.path.join(config["data_dir"], "{sample}_1.fastq.gz"),
    log:
        os.path.join(config["output_dir"], "{sample}", "log", "{sample}_download.log")
    conda:
        DOWNLOAD
    shadow: 
        "minimal" 
    params:
        out_dir = config["data_dir"],
        r2 = os.path.join(config["data_dir"], "{sample}_2.fastq.gz")
    shadow: 
        "minimal"
    shell:
        """
        echo "Starting download for {wildcards.sample}..." > {log}
        # 1. Download data (using current directory via shadow)
        fasterq-dump --split-files --threads {threads} -O . {wildcards.sample} >> {log} 2>&1

        # 2. Compression Loop
        # Check what files were downloaded and compress them
        if [ -f "{wildcards.sample}_1.fastq" ]; then
            gzip "{wildcards.sample}_1.fastq"
            mv "{wildcards.sample}_1.fastq.gz" {output.r1}
        fi

        if [ -f "{wildcards.sample}_2.fastq" ]; then
            gzip "{wildcards.sample}_2.fastq"
            mv "{wildcards.sample}_2.fastq.gz" {params.r2}
        else
            if [ -f "{wildcards.sample}.fastq" ]; then
            gzip "{wildcards.sample}.fastq"
            mv "{wildcards.sample}.fastq.gz" {output.r1}
        fi  
        """