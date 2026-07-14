###################################################################################
#                         workflow/rules/megan.smk                                #
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
#                              version: 06.2026                                   #
###################################################################################

#1. GET MAPPER
rule get_megan6_mapper:
  output:
    db= os.path.join(MEGAN_DIR[0], "megan-map-Feb2022-ue.db"),
  params:
    dir=MEGAN_DIR[0]
  log:
    os.path.join(OUT_DIR, "log", "megan6_download_mapper.log")
  conda:
    DOWNLOAD
  shell:
    """
    wget https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022-ue.db.zip -O {output.db}.zip 2> {log}
    unzip -o {output.db}.zip -d {params.dir} >> {log} 2>&1
    rm {output.db}.zip
    """

#2. READS DIAMOND DAA FORMAT
rule reads_diamond_daa:
  input:
    r1 = get_denovo_r1,
    r2 = get_denovo_r2,
    extra = get_denovo_unpaired
  output:
    hits = temp(os.path.join(OUT_DIR, "{sample}", "megan_reads", "{sample}_reads_report.daa")),
  params:
    db = get_diamond_db_name,
    outfmt = 100,
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity=config["diamond"]["sensitivity"]
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_reads_daa_{sample}.log")
  conda:
    DIAMOND
  shell:
    """
    exec > {log} 2>&1    
    echo "[INFO] Starting Diamond BlastX for Reads..."
    echo "[INFO] Input R1: {input.r1}"
    echo "[INFO] Input R2: {input.r2}"
    echo "[INFO] Input Extra: {input.extra}"
    
    DB_PATH="resources/diamond/{params.db}"
    
    echo "[INFO] Diamond Version:"
    diamond --version
    
    if [ ! -f "$DB_PATH" ] && [ ! -f "$DB_PATH.dmnd" ]; then
        echo "[ERROR] Diamond database not found at: $DB_PATH"
        exit 1
    fi
    echo "[INFO] DB: $DB_PATH"
    
    # Concatenate available inputs into a single temporary file for Diamond
    TMP_READS="{output.hits}.tmp.fastq.gz"
    > "$TMP_READS"
    for f in {input.r1} {input.r2} {input.extra}; do
        if [ -s "$f" ]; then
            cat "$f" >> "$TMP_READS"
        fi
    done
    
     diamond blastx \
      --query "$TMP_READS" \
      --db "$DB_PATH" \
      --out {output.hits} \
      --threads {resources.threads} \
      --outfmt {params.outfmt} \
      --max-target-seqs {params.max_target_seqs} \
      {params.sensitivity}

    rm -f "$TMP_READS"  

    """

rule reads_meganizer:
  input:
    daa = rules.reads_diamond_daa.output.hits,
    db = rules.get_megan6_mapper.output.db
  output:
    daa = temp(os.path.join(OUT_DIR, "{sample}", "megan_reads", "{sample}_reads_report_meganized.daa"))
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "meganizer_{sample}.log")
  conda:
    MEGAN
  shell:
    """
    cp {input.daa} {output.daa}
    
    # Adicionando parâmetros explícitos do LCA para evitar perda de hits
    # -ms 40: Diminui o min-score aceitável (útil para reads curtos/divergentes)
    # -me 0.01: Max e-value aceitável
    # -top 10: Limiar de 10% do melhor score para o LCA
    # -sup 1: Suporte mínimo de 1 read para confirmar um táxon (CRÍTICO)
    
    daa-meganizer \
        -i {output.daa} \
        -mdb {input.db} \
        -ms 40 \
        -me 0.01 \
        -top 10 \
        -sup 1 \
        -t {threads} > {log} 2>&1
    """


rule reads_daa2info:
  input:
    daa = rules.reads_meganizer.output.daa
  output:
    summary = os.path.join(OUT_DIR, "{sample}", "megan_reads", "{sample}_reads_summary.megan")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "daa2info_{sample}.log")
  conda:
    MEGAN
  shell:
    """
    daa2info -i {input.daa} --extractSummaryFile {output.summary} > {log} 2>&1
    """