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
    Adds a clean header.
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
        
        # 2. Filter data        
        grep -v "NOT_FOUND" {input} >> {output.valid_hits} || true
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
        
        # 1. Extrai APENAS os TaxIDs (pulando o cabeçalho)
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # 2. Obtém a linhagem via TaxonKit
        # O output do taxonkit reformat será: TaxID <TAB> Coluna1 <TAB> Coluna2 ...
        taxonkit lineage --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \
        taxonkit reformat --data-dir "${{DB_DIR}}" \
            -f "{{C}}\t{{a}}\t{{d}}\t{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" \
            -F \
             2>> {log} > {output}.lineage.tmp
        
        # 3. Remove a coluna TaxID do output do TaxonKit (pois já temos ela no arquivo original)
        # Ficamos apenas com as colunas de linhagem
        cut -f 2- {output}.lineage.tmp > {output}.reformat_only.tmp
        
        # 4. Verifica integridade (Número de linhas deve bater)
        N_ORIG=$(wc -l < {output}.taxids.tmp)
        N_NEW=$(wc -l < {output}.reformat_only.tmp)
        
        if [ "$N_ORIG" -ne "$N_NEW" ]; then
            echo "[ERROR] Line count mismatch! Input: $N_ORIG, Lineage: $N_NEW" >> {log}
            exit 1
        fi

        # 5. Cria o arquivo final
        # Cabeçalho
        echo -e "{params.base_header}\t{params.lineage_header}" > {output.lineages}
        
        # Cola as colunas originais (1-12) + Taxid (13) + Linhagem
        # Usamos paste para garantir alinhamento tabular perfeito
        tail -n +2 {input.valid_hits} | cut -f 1-13 | \
        paste - {output}.reformat_only.tmp >> {output.lineages}
        
        # Limpeza
        rm {output}.taxids.tmp {output}.lineage.tmp {output}.reformat_only.tmp
        """