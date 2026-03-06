###################################################################################
#                       workflow/rules/diamond.smk                                #
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

rule diamond_blastx_contigs:
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    contigs = get_contigs_path
  output:
    hits = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt"),
    log  = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "diamond.log")
  params:
    db = get_diamond_db_name,
    outfmt = config["diamond"]["outfmt"],
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity = config["diamond"]["sensitivity"]
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_{tool}_{sample}.log")
  conda:
    DIAMOND
  shell:
    """
    exec > {log} 2>&1    
    echo "[INFO] Starting Diamond BlastX for Contigs..."
    echo "[INFO] Input: {input.contigs}"
    
    DB_PATH="resources/diamond/{params.db}"
    
    echo "[INFO] Diamond Version:"
    diamond --version
    
    if [ ! -f "$DB_PATH" ] && [ ! -f "$DB_PATH.dmnd" ]; then
        echo "[ERROR] Diamond database not found at: $DB_PATH"
        exit 1
    fi
    echo "[INFO] DB: $DB_PATH"
    
     diamond blastx \
      --query {input.contigs} \
      --db "$DB_PATH" \
      --out {output.hits} \
      --threads {resources.threads} \
      --outfmt {params.outfmt} \
      --max-target-seqs {params.max_target_seqs} \
      {params.sensitivity} \
      --log
    
    cp {log} {output.log}
    """

rule diamond_blastx_reads:
  input:
    r1 = get_denovo_r1,
    r2 = get_denovo_r2,
    extra = get_denovo_unpaired
  output:
    hits = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_reads_report.txt"),
    log  = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "diamond.log")
  params:
    db = get_diamond_db_name,
    outfmt = config["diamond"]["outfmt"],
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity=config["diamond"]["sensitivity"]
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_reads_{sample}.log")
  conda:
    DIAMOND
  shell:
    """
    # Garante diretório de log aqui também por segurança
    mkdir -p $(dirname {log})

    # 1. Define Temp Dir (Better perfomance)
    LOCAL_TMP=${{SLURM_TMPDIR:-/tmp}}
    MERGED_QUERY="$LOCAL_TMP/{wildcards.sample}_merged_query.fastq.gz"

    echo "[INFO] Using temp dir: $LOCAL_TMP" > {log}

    # 2. Concat in the disk
    echo "[INFO] Concatecating reads..." >> {log}
    cat {input.r1} {input.r2} {input.extra} > "$MERGED_QUERY"

    # 3. Run DIAMOND
    echo "[INFO] Diamond BlastX Reads..." >> {log}
    
    diamond blastx \
      --query "$MERGED_QUERY" \
      --db resources/diamond/{params.db} \
      --out {output.hits} \
      --threads {resources.threads} \
      --outfmt {params.outfmt} \
      --max-target-seqs {params.max_target_seqs} \
      {params.sensitivity} \
      --tmpdir "$LOCAL_TMP" \
      >> {log} 2>&1

    # 4. Limpeza do disco local
    rm -f "$MERGED_QUERY"
    
    # Copia o log final para o output esperado
    cp {log} {output.log}
    """