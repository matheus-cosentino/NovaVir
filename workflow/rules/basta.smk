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
        # Removido directory(), adicionado protected()
        db = protected(os.path.join(BASTA_DB_DIR[0], "prot_mapping.db")),
        gz = protected(os.path.join(BASTA_DB_DIR[0], "prot.accession2taxid.gz"))
    params:
        tax_dir=BASTA_DB_DIR[0]
    conda:
        BASTA
    log:
        os.path.join(OUT_DIR, "log", "basta_download_mapping.log")
    shell:
        """
        # O comando basta download cria os arquivos automaticamente.
        # Se os arquivos output.db e output.gz já existirem, o Snakemake nem entra aqui.
        basta download prot -d {params.tax_dir} > {log} 2>&1
        """

rule basta_download_taxonomy:
    output:
        db    = protected(os.path.join(BASTA_DB_DIR[0], "complete_taxa.db")),
        names = protected(os.path.join(BASTA_DB_DIR[0], "names.dmp")),
        nodes = protected(os.path.join(BASTA_DB_DIR[0], "nodes.dmp"))
    params:
        tax_dir = BASTA_DB_DIR[0]
    conda:
        BASTA
    log:
        os.path.join(OUT_DIR, "log", "basta_download_taxonomy.log")
    shell:
        """
        basta taxonomy -d {params.tax_dir} > {log} 2>&1
        """


rule basta_search:
    input:
        # ancient() evita que pequenas mudanças de data re-disparem essa regra
        mapping_db = ancient(os.path.join(BASTA_DB_DIR[0], "prot_mapping.db")),
        taxonomy   = ancient(os.path.join(BASTA_DB_DIR[0], "complete_taxa.db")),
        query      = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_report.txt")
    output:
        lca=os.path.join(OUT_DIR, "{sample}", "basta_{source}", "{sample}_{source}_lca.tsv"),
        lca_summary=os.path.join(OUT_DIR, "{sample}", "basta_{source}", "{sample}_{source}_lca_summary.tsv")
    params:
        db_type="prot",
        tax_dir=os.path.join(BASTA_DB_DIR[0]),
        algo=config["basta"]["classification"],
        identity=config["basta"]["identity"],
        lengths=config["basta"]["lengths"],
        min_hits= config["basta"]["min_hits"]    
    conda: 
        BASTA
    log:
        os.path.join(OUT_DIR,"{sample}" ,"log", "{sample}_{source}_basta_search.log")
    shell:
        """
        basta sequence {input.query} {output.lca} {params.db_type} \
            -v {output.lca_summary} \
            -d {params.tax_dir} \
            -m {params.min_hits} \
            -i {params.identity} \
            -l {params.lengths} \
            > {log} 2>&1
        """

rule basta_merge_counts:
    input:
        lca_files = get_all_basta_read_outputs
    output:
        table = os.path.join(OUT_DIR, "basta_all", "all_samples_basta_counts.tsv")
    log:
        os.path.join(OUT_DIR, "log", "basta_merge_counts.log")
    run:
        import pandas as pd
        import os
        from collections import defaultdict

        # Dictionary to store counts: counts[Taxonomy][Sample] = Count
        counts = defaultdict(lambda: defaultdict(int))
        all_samples = []

        with open(output.table, 'w') as out_f:
            for lca_file in input.lca_files:
                # Extract sample name
                filename = os.path.basename(lca_file)
                sample_name = filename.replace("_reads_lca.tsv", "")
                all_samples.append(sample_name)
                
                with open(lca_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            # BASTA format: QueryID \t Taxonomy
                            taxonomy = parts[1]
                            # Clean taxonomy string (optional cleanup)
                            taxonomy = taxonomy.replace("_", " ") 
                            counts[taxonomy][sample_name] += 1
            
            # Convert to DataFrame
            df = pd.DataFrame(counts).fillna(0).astype(int)
            
            # The structure is currently Rows=Samples, Cols=Taxa
            # We transpose it to match standard OTU table (Rows=Taxa, Cols=Samples)
            df = df.T
            
            # Save to file
            df.to_csv(output.table, sep='\t', index_label="Taxonomy")

# --- Rule 2: Plot Rarefaction from the Merged Table ---
rule basta_rarefaction_plot:
    input:
        table = os.path.join(OUT_DIR, "basta_all", "all_samples_basta_counts.tsv")
    output:
        pdf = os.path.join(OUT_DIR, "basta_all", "Basta_Rarefaction_Curve.pdf")
    conda:
        R_RAREFACTION
    log:
        os.path.join(OUT_DIR, "log", "basta_rarefaction_plot.log")
    script:
        "../scripts/plot_rarefaction_basta.R"