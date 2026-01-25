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
    Splits the input file into two: one with valid taxids and one with "NOT_FOUND".
    """
    input:
        os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_taxid.tmp")
    output:
        valid_hits = temp(os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_header.tsv")),
        no_lineage = temp(os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_no_lineage.temp"))
    params:
        header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_split_hits.log")
    shell:
        """
        # Create output files with headers
        echo -e "{params.header}" > {output.valid_hits}
        echo -e "{params.header}" > {output.no_lineage}
        
        # CORREÇÃO AQUI: Adicionado BEGIN {{OFS="\\t"}}
        # Isso garante que ao modificar a linha, o awk mantenha os TABs e não troque por espaços.
        
        awk -F'\\t' 'BEGIN {{OFS="\\t"}} 
        {{
            # Limpa carriage return (\r) do Windows se existir
            gsub(/\\r/, "", $NF)
            
            if ($NF == "NOT_FOUND") {{
                print $0 >> "{output.no_lineage}"
            }} else {{
                print $0 >> "{output.valid_hits}"
            }}
        }}' {input}
        """


# 3. Append Lineage to Diamond Results
rule append_lineage:
    input:
        valid_hits = rules.split_hits_by_taxid.output.valid_hits
    output:
        lineages = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_lineage.tsv")
    params:
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tTaxid",
        lineage_header="Lineage\tCelular\tAcelular\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies",
        nodes=os.path.join(BASTA_DB_DIR[0], "nodes.dmp")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_diamond_lineages.log")
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {params.nodes})
        
        # Skip header (criado na regra anterior) e extrai apenas a coluna TaxID (13)
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # 1. Get lineage information
        # TaxonKit retorna: TaxID <TAB> Lineage <TAB> Realm <TAB> Kingdom ...
        taxonkit lineage --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \\
        taxonkit reformat --data-dir "${{DB_DIR}}" \\
            -f "{{C}}\\t{{a}}\\t{{d}}\\t{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" \\
            -F \
             2>> {log} > {output}.lineage.tmp
        
        # 2. Extract reformatted columns only (drop TaxID na col 1 do output do taxonkit)
        cut -f 2- {output}.lineage.tmp > {output}.reformat_only.tmp
        
        # 3. Create final output header
        echo -e "{params.base_header}\\t{params.lineage_header}" > {output.lineages}
        
        # 4. Paste original columns + reformatted lineage
        # Pegamos colunas 1-12 originais + coluna 13 (taxid) original + linhagem nova
        tail -n +2 {input.valid_hits} | cut -f 1-12 | \\
        paste - <(tail -n +2 {input.valid_hits} | cut -f 13) {output}.reformat_only.tmp >> {output.lineages}
        
        # Cleanup
        rm {output}.taxids.tmp {output}.lineage.tmp {output}.reformat_only.tmp
        """