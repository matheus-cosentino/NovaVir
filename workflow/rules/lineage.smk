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
        taxid_map = rules.basta_download_mapping.output.gz
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
      valid_hits=temp(os.path.join(OUT_DIR, "{sample}", "diamond_{source}",  "{sample}_{source}_hits_with_header.tsv"))
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

#################################################
# --- 3. Append Lineage to Diamond Results --- #
################################################ 

rule append_lineage:
    input:
        valid_hits= rules.split_hits_by_taxid.output.valid_hits
    output:
        lineages=os.path.join(OUT_DIR, "{sample}", "diamond_{source}",  "{sample}_{source}_hits_with_lineage.tsv")
    params:
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tTaxid",
        lineage_header="Lineage\tCelular\tAcelular\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies",
        nodes=os.path.join(BASTA_DB_DIR[0], "nodes.dmp"),
        names=os.path.join(BASTA_DB_DIR[0], "names.dmp")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_diamond_lineages.log")
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {params.nodes})
        
        # Skip header and extract taxids
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # 1. Get lineage information, outputting TaxID (default col 1) + 12 reformatted columns (col 2-13), all separated by TAB
        # Usamos '\t' para garantir que os campos sejam separados por tabulação.
        taxonkit lineage --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \\
        taxonkit reformat --data-dir "${{DB_DIR}}" \\
            -f "{{C}}\\t{{a}}\\t{{d}}\\t{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" \\
             2>> {log} > {output}.lineage.tmp
        
        # 2. Extract ONLY the 12 reformatted columns (cols 2-13), dropping the TaxID (col 1),
        # as the TaxID is already available in the main DIAMOND output (col 13).
        # We use cut on the tab-separated intermediate file.
        cut -f 2- {output}.lineage.tmp > {output}.reformat_only.tmp
        
        # 3. Create final output with comprehensive header
        echo -e "{params.base_header}\\t{params.lineage_header}" > {output}
        
        # 4. Combine original DIAMOND data (cols 1-12) + Taxid (col 13) + Reformatted Ranks (cols 14-25)
        tail -n +2 {input.valid_hits} | cut -f 1-12 | \\
        paste - <(tail -n +2 {input.valid_hits} | cut -f 13) {output}.reformat_only.tmp >> {output}
        
        # Cleanup
        rm {output}.taxids.tmp {output}.lineage.tmp {output}.reformat_only.tmp
        """