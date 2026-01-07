###################################################################################
#                       workflow/rules/lineage.smk                                #
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


#################################################
# --- 1. Map Proteins Id Diamond for Taxid --- #
################################################ 

rule map_accession_to_taxid:
    """
    Maps protein IDs (Subject ID, column 2 of DIAMOND) to TaxIDs.
    """
    input:
        hit_file = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_report.txt"),
        taxid_map = rules.basta_download_mapping.output
    output:
        ids = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_taxid.tmp")
    shadow: 
        "minimal"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_map_acc_prot.log")
    script:
        "../scripts/map_acc_to_taxid.py"


##################################################
# --- 2. Split Hits for Taxid and Not Found --- #
################################################# 

rule split_hits_by_taxid:
    """
    Filters the input file to keep only hits with valid TaxIDs.
    """
    input:
      os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_taxid.tmp")
    output:
      valid_hits=os.path.join(OUT_DIR, "{sample}", "diamond_{source}",  "{sample}_{source}_hits_with_lineage.tsv")
    params:
      header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    log:
      os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_split_hits.log")
    shell:
        """
        # 1. Create the output file with the header
        echo -e "{params.header}" > {output.valid_hits}
        
        # 2. Filter data: Skip header (NR==1) AND only print if last column is NOT "NOT_FOUND"
        awk -F'\\t' '
        NR==1 {{ next }} 
        $NF != "NOT_FOUND" {{
            print $0 >> "{output.valid_hits}"
        }}' {input}
        """


