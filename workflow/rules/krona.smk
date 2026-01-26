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
        acc_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted")
    input:
        names = os.path.join(BASTA_DB_DIR[0], "names.dmp"),
        nodes = os.path.join(BASTA_DB_DIR[0], "nodes.dmp"),
        acc2tax = os.path.join(BASTA_DB_DIR[0], "prot.accession2taxid.gz")
    params:
        tax_dir=KRONA_DB_DIR[0]
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "log", "krona_update_taxonomy.log")
    shell:
        """
        echo "[INFO] Configurando diretórios..." > {log}
        TARGET_DIR=$(readlink -f {params.tax_dir})
        
        # 1. Taxonomia (Isso funciona, vamos manter)
        echo "[INFO] Linkando arquivos de taxonomia..." >> {log}
        ln -sf $(readlink -f {input.names}) $TARGET_DIR/names.dmp
        ln -sf $(readlink -f {input.nodes}) $TARGET_DIR/nodes.dmp
        
        echo "[INFO] Executando ktUpdateTaxonomy.sh..." >> {log}
        ktUpdateTaxonomy.sh --only-build $TARGET_DIR >> {log} 2>&1

        # 2. Accessions (Fazer Manualmente = 100% Seguro)
        # O script original falha em achar o arquivo, então nós mesmos criamos o output.
        # Pegamos colunas 2 (Accession.Version) e 3 (TaxID) e ordenamos.
        
        echo "[INFO] Gerando accession2taxid.sorted manualmente..." >> {log}
        
        zcat {input.acc2tax} | \
        cut -f 2,3 | \
        sort -k 1,1 > {output.acc_sorted} 2>> {log}
        
        echo "[INFO] Concluído." >> {log}
        """

rule krona_kraken2:
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_output.txt"),
        tax_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_{tool}_kraken2_krona.html")
    params:
        tax_dir=KRONA_DB_DIR[0]
    conda:
        KRONA
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_kraken2_{tool}_{sample}.log")
    shell:
        """
        ktImportTaxonomy \
            -tax {params.tax_dir} \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
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
        TMP_DB="{wildcards.sample}_{wildcards.read_type}_reads_krona_db"
        mkdir -p $TMP_DB
        ln -sf $(readlink -f {input.tax_sorted}) $TMP_DB/accession2taxid.sorted
        ln -sf $(readlink -f {input.tax_tab}) $TMP_DB/taxonomy.tab

        ktImportTaxonomy \
            -tax $TMP_DB \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
        
        rm -rf $TMP_DB
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
        TMP_DB="{wildcards.sample}_reads_diamond_krona_db"
        mkdir -p $TMP_DB
        ln -sf $(readlink -f {input.tax_sorted}) $TMP_DB/accession2taxid.sorted
        ln -sf $(readlink -f {input.tax_tab}) $TMP_DB/taxonomy.tab

        ktImportBLAST \
            {input.report} \
            -tax $TMP_DB \
            -o {output.html} \
            > {log} 2>&1
            
        rm -rf $TMP_DB
        """

rule krona_diamond_contigs:
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt"),
        # Estes arquivos já foram gerados e estão na raiz do KRONA_DB_DIR
        tax_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_{tool}_diamond_krona.html")
    conda:
        KRONA
    params:
        tax_dir=KRONA_DB_DIR[0]
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "krona_diamond_contigs_{tool}_{sample}.log")
    shell:
        """
        echo "[INFO] Iniciando setup do Krona..." > {log}

        # 1. Definir caminhos absolutos
        DB_ROOT=$(readlink -f {params.tax_dir})
        ACC_SUBDIR="$DB_ROOT/accession2taxid"
        SRC_FILE="$DB_ROOT/accession2taxid.sorted"

        # 2. Validar existência dos arquivos originais
        if [ ! -f "$DB_ROOT/taxonomy.tab" ]; then
            echo "[ERROR] taxonomy.tab não encontrado em $DB_ROOT" >> {log}
            exit 1
        fi
        
        # 3. Criar diretório de accessions (se não existir)
        mkdir -p "$ACC_SUBDIR"

        echo "[INFO] Criando links universais em $ACC_SUBDIR..." >> {log}

        # 4. A ESTRATÉGIA DE ADAPTAÇÃO:
        # O Krona pode procurar por qualquer um desses nomes. 
        # Nós criamos links de TODOS eles apontando para o seu arquivo único gerado.
        
        FILES_TO_LINK=(
            "accession2taxid.sorted"
            "prot.accession2taxid.sorted"
            "nucl_gb.accession2taxid.sorted"
            "nucl_wgs.accession2taxid.sorted"
            "nucl_est.accession2taxid.sorted"
            "nucl_gss.accession2taxid.sorted"
            "dead_prot.accession2taxid.sorted"
            "dead_nucl.accession2taxid.sorted"
            "dead_wgs.accession2taxid.sorted"
        )

        for FILE_NAME in "${{FILES_TO_LINK[@]}}"; do
            # Remove link antigo se existir e cria um novo (hard link preferencialmente)
            ln -f "$SRC_FILE" "$ACC_SUBDIR/$FILE_NAME" 2>/dev/null || ln -sf "$SRC_FILE" "$ACC_SUBDIR/$FILE_NAME"
        done

        # Debug: Mostrar o que foi criado
        ls -l "$ACC_SUBDIR" >> {log}

        # 5. Executar ktImportBLAST
        echo "[INFO] Executando ktImportBLAST..." >> {log}
        
        # LC_ALL=C é CRÍTICO para o Perl ler arquivos sorted corretamente
        export LC_ALL=C
        
        ktImportBLAST \
            {input.report} \
            -tax "$DB_ROOT" \
            -o {output.html} \
            >> {log} 2>&1
        """