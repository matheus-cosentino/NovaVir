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
        db_dir = directory(KRONA_DB_DIR[0]),
        marker = touch(os.path.join(KRONA_DB_DIR[0], "accessions.done"))
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
        mkdir -p {output.db_dir}
        echo "[INFO] Linking taxonomy files..." > {log}
        ln -sf $(readlink -f {input.names}) {output.db_dir}/names.dmp
        ln -sf $(readlink -f {input.nodes}) {output.db_dir}/nodes.dmp

        echo "[INFO] Building Taxonomy Tree..." >> {log}
        ktUpdateTaxonomy.sh --only-build {output.db_dir} >> {log} 2>&1

        echo "[INFO] Processing Accessions..." >> {log}
        ln -sf $(readlink -f {input.acc2tax}) {output.db_dir}/accession2taxid.gz

        echo "[INFO] Finding Krona updateAccessions..." >> {log}
        KRONA_BIN_DIR=$(dirname $(which ktImportBLAST))
        HIDDEN_SCRIPT="$KRONA_BIN_DIR/../opt/krona/updateAccessions.sh"

        if command -v ktUpdateAccessions.sh &> /dev/null; then
            ktUpdateAccessions.sh --file {input.acc2tax} --taxonomy {output.db_dir} >> {log} 2>&1
        elif [ -f "$HIDDEN_SCRIPT" ]; then
            echo "[INFO] Found hidden script at $HIDDEN_SCRIPT" >> {log}
            bash "$HIDDEN_SCRIPT" --file {input.acc2tax} --taxonomy {output.db_dir} >> {log} 2>&1
        else
            echo "[WARN] Could not find updateAccessions.sh anywhere." >> {log}
        fi
        """

rule krona_kraken2:
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_output.txt"),
        # IMPORTANTE: Adiciona o DB como input para o Snakemake montar no shadow dir
        tax_db = KRONA_DB_DIR[0]
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_{tool}_kraken2_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_kraken2_{tool}_{sample}.log")
    shell:
        """
        ktImportTaxonomy \
            -tax {input.tax_db} \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
        """

rule krona_reads_kraken:
    wildcard_constraints:
        read_type="paired|unpaired"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_{read_type}_reads_output.txt"),
        tax_db = KRONA_DB_DIR[0]
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_reads", "{sample}_{read_type}_kraken2_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_kraken_reads_{read_type}_{sample}.log")
    shell:
        """
        ktImportTaxonomy \
            -tax {input.tax_db} \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
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
        tax_db = KRONA_DB_DIR[0]
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_reads", "{sample}_reads_diamond_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_diamond_reads_{sample}.log")
    shell:
        """
        ktImportBLAST \
            {input.report} \
            -tax {input.tax_db} \
            -o {output.html} \
            > {log} 2>&1
        """

rule krona_diamond_contigs:
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt"),
        # CORREÇÃO PRINCIPAL: Banco de dados adicionado como INPUT
        tax_db = KRONA_DB_DIR[0]
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_{tool}_diamond_krona.html")
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_diamond_contigs_{tool}_{sample}.log")
    shell:
        """
        ktImportBLAST \
            {input.report} \
            -tax {input.tax_db} \
            -o {output.html} \
            > {log} 2>&1
        """