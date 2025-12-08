###################################################################################
#                       workflow/rules/kraken2.smk                                #
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

rule kraken2_contigs:
    input:
      contigs= get_contigs_path
    output:
      report= "{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_report.txt"
      #out="{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_output.txt"
    params:
      #db=f"{workflow.basedir}/{config["resources"]["kraken2"]}",
      db=config["resources"]["kraken2"],        
      confidence=config["kraken2"]["confidence"]
    log:
      "{out_dir}/{sample}/log/kraken2_contigs_{tool}_{sample}.log"
    conda:
      KRAKEN2        
    shell:
      """
      kraken2 \
      --db {params.db} \
      --confidence {params.confidence} \
      --report {output.report} \
      --threads {resources.threads} \
      --memory-mapping \
      {input.contigs} \
      > {log} 2>&1
      """

rule kraken_biom_contig:
    input:
      report= "{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_report.txt",
    output:
      biom= "{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_biom.txt",
    log:
      "{out_dir}/{sample}/log/kraken2_contigs_{tool}_{sample}_biom.log"
    params:
      maximun=config["kraken_biom"]["max"],
      minimum=config["kraken_biom"]["min"],
      out_format=config["kraken_biom"]["format"]
    conda:
        KRAKEN2_BIOM
    shell:
        """
        kraken-biom \
        {input.report} \
        --max {params.maximun} \
        --min {params.minimum} \
        --fmt {params.out_format} \
        -o {output.biom} \
        > {log} 2>&1
        """

rule kraken2_reads_paired:
    input:
     r1 = os.path.join(OUT_DIR, "{sample}", "trimmed", "{sample}_1.fastq.gz"),
     r2 = os.path.join(OUT_DIR, "{sample}", "trimmed" , "{sample}_2.fastq.gz"),
     extra = os.path.join(OUT_DIR, "{sample}", "trimmed", "{sample}_orphans.fastq.gz")
    output:
      report= "{out_dir}/{sample}/kraken2_reads/{sample}_paired_reads_report.txt"
      #out="{out_dir}/{sample}/kraken2_reads/{sample}_paired_reads_output.txt"
    params:
      #db=f"{workflow.basedir}/{config["resources"]["kraken2"]}",
      db=config["resources"]["kraken2"],        
      confidence=config["kraken2"]["confidence"]
    log:
      "{out_dir}/{sample}/log/kraken2_paired_reads_{sample}.log"
    conda:
      KRAKEN2        
    shell:
      """
      kraken2 \
      --db {params.db} \
      --confidence {params.confidence} \
      --report {output.report} \
      --memory-mapping \
      --threads {resources.threads} \
      {input.r1} {input.r2} {input.extra} \
      > {log} 2>&1
      """

rule kraken2_reads_unpaired:
    input:   
      r1= os.path.join(OUT_DIR, "{sample}", "trimmed", "{sample}_unp.fastq.gz")
    output:
      report= "{out_dir}/{sample}/kraken2_reads/{sample}_unpaired_reads_report.txt"
      #out="{out_dir}/{sample}/kraken2_reads/{sample}_unpared_reads_output.txt"
    params:
      #db=f"{workflow.basedir}/{config["resources"]["kraken2"]}",
      db=config["resources"]["kraken2"],        
      confidence=config["kraken2"]["confidence"]
    log:
      "{out_dir}/{sample}/log/kraken2_paired_reads_{sample}.log"
    conda:
      KRAKEN2        
    shell:
      """
      kraken2 \
      --db {params.db} \
      --confidence {params.confidence} \
      --report {output.report} \
      --output {output.out} \
      --memory-mapping \
      --threads {resources.threads} \
      {input.r1} \
      > {log} 2>&1
      """

rule kraken_biom_reads:
    input:
      report= "{out_dir}/{sample}/kraken2_reads/{sample}_{paired}_reads_report.txt"
    output:
      biom= "{out_dir}/{sample}/kraken2_reads/{sample}_{paired}_reads_biom.txt",
    params:
      maximun=config["kraken_biom"]["max"],
      minimum=config["kraken_biom"]["min"],
      out_format=config["kraken_biom"]["format"]
    log:
        "{out_dir}/{sample}/log/kraken2_biom_{paired}_reads_{sample}.log"
    conda:
        KRAKEN2_BIOM
    shell:
        """
        kraken-biom \
        {input.report} \
        --max {params.maximun} \
        --min {params.minimum} \
        --fmt {params.out_format} \
        -o {output.biom} \
        > {log} 2>&1
        """