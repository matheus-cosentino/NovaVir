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

rule basta_prepare_mapping:
    input:
        mapping_file = config["resources"]["taxonmap"]
    output:
        temp(os.path.join(OUT_DIR, "temp", "basta_mapping.txt"))
    log:
        os.path.join(OUT_DIR, "log", "basta_prepare_mapping.log")
    shell:
        """
        if [[ "{input.mapping_file}" == *.gz ]]; then
            gunzip -c {input.mapping_file} > {output} 2> {log}
        else
            ln -sfr {input.mapping_file} {output} 2> {log}
        fi
        """

rule basta_createdb:
    input:
        mapping=os.path.join(OUT_DIR, "temp", "basta_mapping.txt")
    output:
        directory(os.path.join("resources", "basta_db"))
    params:
        acc_col=1,
        taxid_col=2,
        db_name="prot_mapping"
    conda:
        BASTA
    log:
        os.path.join(OUT_DIR, "log", "basta_createdb.log")
    shell:
        """
        # Create the DB with a specific name
        basta create_db {input.mapping} {params.db_name} {params.acc_col} {params.taxid_col} -d {output} >> {log} 2>&1
        """

rule basta_search:
    input:
        mapping_db=directory(os.path.join("resources", "basta_db")),
        query=os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_report.txt"),
      
    output:
        lca=os.path.join(OUT_DIR, "{sample}", "basta_{source}", "{sample}_{source}_lca.tsv")
    params:
        db_type="prot_mapping",
        tax_dir=os.path.dirname(config["resources"]["taxonnodes"][0])
    conda: BASTA
    log:
        os.path.join(OUT_DIR,"{sample}" ,"log", "{sample}_{source}_basta_search.log")
    shell:
        """
        basta sequence {input.query} {output.lca} {params.db_type} \
            -d {input.mapping_db} \
            > {log} 2>&1
        """
