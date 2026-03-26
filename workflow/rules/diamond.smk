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
    db_name = get_diamond_db_name,
    db = lambda wildcards: os.path.join(DIAMOND_DIR[0], get_diamond_db_name(wildcards)),
    outfmt = config["diamond"]["outfmt"],
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity = config["diamond"]["sensitivity"]
  conda:
    DIAMOND
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_{tool}_{sample}.log")
  shell:
    """
    exec > {log} 2>&1    
    echo "[INFO] Starting Diamond BlastX for Contigs..."
    echo "[INFO] Input: {input.contigs}"
    
    DB_PATH="{params.db}"
    
    echo "[INFO] Diamond Version:"
    diamond --version
    
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
    db_name = get_diamond_db_name,
    db = lambda wildcards: os.path.join(DIAMOND_DIR[0], get_diamond_db_name(wildcards)),
    outfmt = config["diamond"]["outfmt"],
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity=config["diamond"]["sensitivity"]
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_reads_{sample}.log")
  conda: 
    DIAMOND
  shell:
    """
    exec > {log} 2>&1    
    echo "[INFO] Starting Diamond BlastX for Reads..."
    echo "[INFO] Input R1: {input.r1}"
    echo "[INFO] Input R2: {input.r2}"
    echo "[INFO] Input Extra: {input.extra}"
    
    
    
    DB_PATH="{params.db}"
    echo "[INFO] DB: $DB_PATH"
    
    echo "[INFO] Diamond Version:"
    diamond  --version
        
    # Concatenate available inputs into a single temporary file for Diamond
    TMP_READS="{output.hits}.tmp.fastq.gz"
    > "$TMP_READS"
    for f in {input.r1} {input.r2} {input.extra}; do
        if [ -s "$f" ]; then
            cat "$f" >> "$TMP_READS"
        fi
    done
    
     diamond  blastx \
      --query "$TMP_READS" \
      --db "$DB_PATH" \
      --out {output.hits} \
      --threads {resources.threads} \
      --outfmt {params.outfmt} \
      --max-target-seqs {params.max_target_seqs} \
      {params.sensitivity} \
      --log
    
    cp {log} {output.log}
    rm -f "$TMP_READS"
    
    """