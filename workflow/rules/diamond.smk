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

###################################################################################
#                       workflow/rules/diamond.smk                                #
#                         MSc. Matheus Cosentino                                  # 
###################################################################################
#                              version: 12.2025                                   #
###################################################################################

rule diamond_blastx_contigs:
  #wildcard_constraints:
    #tool="^(?!reads$).*$"
  input:
    contigs = get_contigs_path
  output:
    hits = os.path.join(OUT_DIR, "{sample}", "diamond_contigs", "{sample}_contigs_report.txt"),
    log  = os.path.join(OUT_DIR, "{sample}", "diamond_contigs", "diamond.log")
  params:
    #tool = config["tool"]["denovo"],
    db = get_diamond_db_name,
    outfmt = config["diamond"]["outfmt"],
    max_target_seqs = config["diamond"]["max_target_seqs"],
    evalue = config["diamond"]["evalue"]
  log:
    # This log uses {tool}, so output MUST also use {tool}
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_contigs_{sample}.log")
  conda:
    DIAMOND
  shell:
    """
    diamond blastx \
    --query {input.contigs} \
    --db resources/diamond/{params.db} \
    --out {output.hits} \
    --threads {resources.threads} \
    --outfmt {params.outfmt} \
    --max-target-seqs {params.max_target_seqs} \
    --evalue {params.evalue} \
    --log \
    &> {log}

    mv diamond.log {output.log}
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
    evalue = config["diamond"]["evalue"]
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "diamond_reads_{sample}.log")
  conda:
    DIAMOND
  shell:
    """
    diamond blastx \
      --query {input.r1} {input.r2} {input.extra} \
      --db resources/diamond/{params.db} \
      --out {output.hits} \
      --threads {resources.threads} \
      --outfmt {params.outfmt} \
      --max-target-seqs {params.max_target_seqs} \
      --evalue {params.evalue} \
      --log \
      &> {log}
          
    mv diamond.log {output.log}
    """