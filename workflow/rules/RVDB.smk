###################################################################################
#                        workflow/rules/RVDB.smk                                  #
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
#                              version: 03.2026                                   #
###################################################################################

rule get_rvdb:
  output:
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.xz"),
    sql =  os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite.xz")
  log:
    os.path.join(OUT_DIR, "log", "get_rvdb.log")
  conda:
    DOWNLOAD
  shell:
    """
    curl -L https://rvdb-prot.pasteur.fr/files/U-RVDBv31.0-prot.hmm.xz -o  {output.hmm} 2>> {log}
    curl -L https://rvdb-prot.pasteur.fr/files/U-RVDBv31.0-prot-hmm.sqlite.xz -o  {output.sql} 2>> {log}
    """

rule prepare_rvdb:
  input:
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.xz"),
    sql = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite.xz")
  output:
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm"),
    sql = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite"),
    i1 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3m"),
    i2 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3i"),
    i3 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3f"),
    i4 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3p")
  log:
    os.path.join(OUT_DIR, "log", "prepare_rvdb.log")
  conda:
    HMMER
  shell:
    """
    xz -dc {input.hmm} > {output.hmm} 2>> {log}
    xz -dc {input.sql} > {output.sql} 2>> {log}
    hmmpress -f {output.hmm} >> {log} 2>&1
    """

rule rvdb_search:
  input:
    orfs = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta"),
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm"),
    indices = multiext(os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm"), ".h3m", ".h3i", ".h3f", ".h3p")
  output:
    result = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_RVDB_results.tsv")
  params:
    evalue=config['hmmscan']['evalue']
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "rvdb_search_{tool}_{sample}.log")
  conda:
    HMMER
  shell:
    """
    hmmscan --cpu {resources.threads} -E {params.evalue} --tblout {output.result} {input.hmm} {input.orfs} 2>> {log}
    """   

rule rvdb_summarize:
  input:
    tbl = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_RVDB_results.tsv"),
    sqlite = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite")
  output:
    csv = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_RVDB_Summary.csv")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "rvdb_summarize_{tool}_{sample}.log")
  conda:
    CORE
  script:
    "../scripts/rvdb_summarize.py"