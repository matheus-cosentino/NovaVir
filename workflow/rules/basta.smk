###################################################################################
#                         workflow/rules/basta.smk                                # 
#                           MSc. Matheus Cosentino                                # 
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

rule basta_download_mapping:
    output:
        db = directory(os.path.join(BASTA_DB_DIR[0], "prot_mapping.db")),
        gz=os.path.join(BASTA_DB_DIR[0], "prot.accession2taxid.gz")
    params:
        tax_dir=os.path.join(BASTA_DB_DIR[0])
    conda:
        BASTA
    log:
        os.path.join(OUT_DIR, "log", "basta_download_mapping.log")
    shell:
        "basta download prot -d {params.tax_dir} > {log} 2>&1"

rule basta_download_taxonomy:
    output:
        gz=directory(os.path.join(BASTA_DB_DIR[0], "complete_taxa.db"))
    params:
        tax_dir=os.path.join(BASTA_DB_DIR[0])
    conda:
        BASTA
    log:
        os.path.join(OUT_DIR, "log", "basta_download_taxonomy.log")
    shell:
        "basta taxonomy -d {params.tax_dir} > {log} 2>&1"


rule basta_search:
    input:
        mapping_db=os.path.join(BASTA_DB_DIR[0], "prot_mapping.db"),
        query=os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_report.txt"),
        taxonomy=os.path.join(BASTA_DB_DIR[0], "complete_taxa.db")
    output:
        lca=os.path.join(OUT_DIR, "{sample}", "basta_{source}", "{sample}_{source}_lca.tsv"),
        lca_summary=os.path.join(OUT_DIR, "{sample}", "basta_{source}", "{sample}_{source}_lca_summary.tsv")
    params:
        db_type="prot",
        tax_dir=os.path.join(BASTA_DB_DIR[0]),
        algo=config["basta"]["classification"]
    conda: 
        BASTA
    log:
        os.path.join(OUT_DIR,"{sample}" ,"log", "{sample}_{source}_basta_search.log")
    shell:
        """
        basta sequence {input.query} {output.lca} {params.db_type} \
            -v {output.lca_summary} \
            -d {params.tax_dir} \
            -m 1 \
            > {log} 2>&1
        """
