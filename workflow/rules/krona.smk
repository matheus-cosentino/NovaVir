###################################################################################
#                         workflow/rules/krona.smk                                #
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
#                              version: 01.2026                                   #
###################################################################################

rule krona_update_taxonomy:
    output:
        tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab"),
        sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted")
    input:
        names = os.path.join(BASTA_DB_DIR[0], "names.dmp"),
        nodes = os.path.join(BASTA_DB_DIR[0], "nodes.dmp"),
        acc2tax = os.path.join(BASTA_DB_DIR[0], "prot.accession2taxid.gz")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "log", "krona_update_taxonomy.log")
    shell:
        """
        mkdir -p $(dirname {output.tab})

        echo "[INFO] Linking taxonomy files..." > {log}
        ln -sf $(readlink -f {input.names}) $(dirname {output.tab})/names.dmp
        ln -sf $(readlink -f {input.nodes}) $(dirname {output.tab})/nodes.dmp

        echo "[INFO] Building Taxonomy Tree..." >> {log}
        ktUpdateTaxonomy.sh --only-build $(dirname {output.tab}) >> {log} 2>&1

        echo "[INFO] Generating accession2taxid.sorted manually..." >> {log}
        # Criação manual do arquivo sorted (mais robusto que o script interno do Krona)
        zcat {input.acc2tax} | \
        sed '1d' | \
        cut -f 2,3 | \
        LC_ALL=C sort -k1,1 --parallel={threads} -S 10G > {output.sorted} 2>> {log}
        
        if [ ! -s {output.sorted} ]; then
             echo "[ERROR] Failed to create sorted file." >> {log}
             exit 1
        fi
        
        echo "[INFO] Done." >> {log}
        """

rule krona_kraken2:
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_output.txt"),
        # Input explícito para garantir que o Snakemake monte os arquivos no shadow dir
        tax_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_{tool}_kraken2_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_kraken2_{tool}_{sample}.log")
    shell:
        """
        # Estratégia de DB Local Temporário
        mkdir -p tmp_krona_db
        ln -sf $(readlink -f {input.tax_sorted}) tmp_krona_db/accession2taxid.sorted
        ln -sf $(readlink -f {input.tax_tab}) tmp_krona_db/taxonomy.tab

        ktImportTaxonomy \
            -tax tmp_krona_db \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
            
        rm -rf tmp_krona_db
        """

rule krona_reads_kraken:
    wildcard_constraints:
        read_type="paired|unpaired"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_{read_type}_reads_output.txt"),
        tax_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_reads", "{sample}_{read_type}_kraken2_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_kraken_reads_{read_type}_{sample}.log")
    shell:
        """
        mkdir -p tmp_krona_db
        ln -sf $(readlink -f {input.tax_sorted}) tmp_krona_db/accession2taxid.sorted
        ln -sf $(readlink -f {input.tax_tab}) tmp_krona_db/taxonomy.tab

        ktImportTaxonomy \
            -tax tmp_krona_db \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
        
        rm -rf tmp_krona_db
        """

rule krona_basta:
    input:
        lca = os.path.join(OUT_DIR, "{sample}", "basta_{source}", "{sample}_{source}_lca.tsv")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{source}", "{sample}_{source}_basta_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_basta_{source}_{sample}.log")
    shell:
        """
        cut -f 2 {input.lca} | \
        sed 's/; /\\t/g' | \
        ktImportText -o {output.html} - \
        > {log} 2>&1
        """

rule krona_diamond_reads:
    input:
        report = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_reads_report.txt"),
        tax_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_reads", "{sample}_reads_diamond_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_diamond_reads_{sample}.log")
    shell:
        """
        mkdir -p tmp_krona_db
        ln -sf $(readlink -f {input.tax_sorted}) tmp_krona_db/accession2taxid.sorted
        ln -sf $(readlink -f {input.tax_tab}) tmp_krona_db/taxonomy.tab

        ktImportBLAST \
            {input.report} \
            -tax tmp_krona_db \
            -o {output.html} \
            > {log} 2>&1
            
        rm -rf tmp_krona_db
        """

rule krona_diamond_contigs:
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt"),
        tax_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_{tool}_diamond_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_diamond_contigs_{tool}_{sample}.log")
    shell:
        """
        # CRÍTICO: Criar DB local linkado para enganar o Krona
        # Isso resolve problemas de caminhos absolutos vs shadow directories
        mkdir -p tmp_krona_db
        ln -sf $(readlink -f {input.tax_sorted}) tmp_krona_db/accession2taxid.sorted
        ln -sf $(readlink -f {input.tax_tab}) tmp_krona_db/taxonomy.tab

        ktImportBLAST \
            {input.report} \
            -tax tmp_krona_db \
            -o {output.html} \
            > {log} 2>&1
        
        rm -rf tmp_krona_db
        """