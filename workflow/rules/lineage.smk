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

# Get rarefaction taxonomic level, with a default value
RAREFACTION_TAX_LEVEL = config.get('rarefaction', {}).get('tax_level', 'Genus')

#################################################
# --- 1. Map Proteins Id Diamond for Taxid --- #
################################################ 
rule download_prot:
    """
    Download Protein acc to taxid data
    """
    output:
        gz = protected(os.path.join(TAXONOMY_DIR[0], "prot.accession2taxid.gz")),
        nodes = protected(os.path.join(TAXONOMY_DIR[0], "nodes.dmp")),
        names = protected(os.path.join(TAXONOMY_DIR[0], "names.dmp"))
    log:
        os.path.join(OUT_DIR, "log", "tax_download_mapping.log")
    params:
        taxon_gz = protected(os.path.join(TAXONOMY_DIR[0], "taxdump.tar.gz")),
        #tax_dir = TAXONOMY_DIR[0]
    conda:
        DOWNLOAD
    shell:
        """
        wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz -O {output.gz} 2> {log}
        wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O {params.taxon_gz} 2>> {log}
        tar -zxvf {params.taxon_gz} -C {TAXONOMY_DIR[0]} names.dmp nodes.dmp 2>> {log}
        """


#################################################
# --- 2. Map Proteins Id Diamond for Taxid --- #
################################################ 

rule map_accession_to_taxid:
    """
    Maps protein IDs (Subject ID, column 2 of DIAMOND) to TaxIDs.
    """
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled|reads"
    input:
        hit_file = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt"),
        taxid_map = ancient(rules.download_prot.output.gz)
    output:
        ids = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_taxid.tmp")
    shadow: 
        "minimal"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{tool}_map_acc_prot.log")
    script:
        "../scripts/map_acc_to_taxid.py"


##################################################
# --- 3. Split Hits for Taxid and Not Found --- #
################################################# 
rule split_hits_by_taxid:
    """
    Splits the input file into two: one with valid taxids and one with "NOT_FOUND".
    """
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled|reads"
    input:
        os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_taxid.tmp")
    output:
        valid_hits = temp(os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_header.tsv")),
        no_lineage = temp(os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_no_lineage.temp"))
    params:
        header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{tool}_split_hits.log")
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


##################################################
# --- 3. Split Hits for Taxid and Not Found --- #
################################################# 
rule append_lineage:
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled|reads"
    input:
        valid_hits = rules.split_hits_by_taxid.output.valid_hits,
        nodes = ancient(os.path.join(TAXONOMY_DIR[0], "nodes.dmp"))
    output:
        lineages = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_lineage.tsv")
    params:
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tTaxid",
        lineage_header="Lineage\tCelular\tAcelular\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{tool}_diamond_lineages.log")
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {input.nodes})
        
        # Skip header (criado na regra anterior) e extrai apenas a coluna TaxID (13)
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # 1. Get lineage information
        # TaxonKit retorna: TaxID <TAB> Lineage <TAB> Realm <TAB> Kingdom ...
        taxonkit lineage --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \
        taxonkit reformat --data-dir "${{DB_DIR}}" \
            -f "{{C}}\\t{{a}}\\t{{d}}\\t{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" \
            -F \
             2>> {log} > {output}.lineage.tmp
        
        # 2. Extract reformatted columns only (drop TaxID na col 1 do output do taxonkit)
        cut -f 2- {output}.lineage.tmp > {output}.reformat_only.tmp
        
        # 3. Create final output header
        echo -e "{params.base_header}\\t{params.lineage_header}" > {output.lineages}
        
        # 4. Paste original columns + reformatted lineage
        # Pegamos colunas 1-12 originais + coluna 13 (taxid) original + linhagem nova
        tail -n +2 {input.valid_hits} | cut -f 1-12 | \
        paste - <(tail -n +2 {input.valid_hits} | cut -f 13) {output}.reformat_only.tmp >> {output.lineages}
        
        # Cleanup
        rm {output}.taxids.tmp {output}.lineage.tmp {output}.reformat_only.tmp
        """

# --- Rule to merge diamond counts from reads ---
rule diamond_merge_counts:
    input:
        lineage_files = get_all_diamond_lineage_for_reads
    output:
        table = os.path.join(OUT_DIR, "diamond_all", f"all_samples_diamond_counts_by_{RAREFACTION_TAX_LEVEL}.tsv")
    params:
        tax_level = RAREFACTION_TAX_LEVEL
    log:
        os.path.join(OUT_DIR, "log", "diamond_merge_counts.log")
    run:
        import pandas as pd
        import os
        from collections import defaultdict

        tax_level = params.tax_level
        counts = defaultdict(lambda: defaultdict(int))
        
        with open(output.table, 'w') as out_f:
            for lineage_file in input.lineage_files:
                sample_name = os.path.basename(lineage_file).replace("_reads_hits_with_lineage.tsv", "")
                
                df_hits = pd.read_csv(lineage_file, sep='\t')
                
                if tax_level not in df_hits.columns:
                    raise ValueError(f"Taxonomic level '{tax_level}' not found in columns of {lineage_file}")

                tax_counts = df_hits[tax_level].value_counts()
                
                for tax_name, count in tax_counts.items():
                    if pd.notna(tax_name) and tax_name != '-':
                        counts[tax_name][sample_name] += count
            
            df_counts = pd.DataFrame(counts).fillna(0).astype(int)
            
            if not df_counts.empty:
                df_counts = df_counts.T
            
            df_counts.to_csv(output.table, sep='\t', index_label="Taxonomy")

# --- Rule to plot rarefaction from the merged Diamond table ---
rule diamond_rarefaction_plot:
    input:
        table = rules.diamond_merge_counts.output.table
    output:
        pdf = os.path.join(OUT_DIR, "diamond_all", "Diamond_Rarefaction_Curve.pdf")
    params:
        title = f"Diamond Rarefaction Curves (by {RAREFACTION_TAX_LEVEL})"
    conda:
        R_RAREFACTION
    log:
        os.path.join(OUT_DIR, "log", "diamond_rarefaction_plot.log")
    script:
        "../scripts/plot_rarefaction.R"
