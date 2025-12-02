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
      report= "{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_report.txt",
      out="{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_output.txt"
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
      --report {output.hits} \
      --output {output.out} \
      --threads {resources.threads} \
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
    conda:
        KRAKEN2
    shell:
        """
        kraken-biom \
        {input.report}
        --min F \
        --fmt tsv \
        -o {output.biom} \
        > {log} 2>&1
        """