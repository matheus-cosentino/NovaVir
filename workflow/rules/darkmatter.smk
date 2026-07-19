###################################################################################
#                      workflow/rules/duskmatter.smk                              #
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


######################################################
# --- 1. Filter Contig Fasta With no Diamond Hit --- #
###################################################### 
rule get_nohit_fasta:
  """
    Get contigs within no hits against the nr Diamond database
  """
  # ADICIONAR ESTA CONSTRAINT
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    fasta = get_contigs_path,
    diamond = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_lineage.tsv")  
  output:
    nohits = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_{tool}_nohit.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "darkmatter_{tool}_{sample}_nohits.log")
  conda:
    CORE
  script:
    "../scripts/filter_nohits.py"


#################################################
# --- 2. Find Orfs for the no Hits Contigs --- #
################################################ 
rule find_orfs:
  input: 
    fasta = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_{tool}_nohit.fasta")
  output:
    orfs = temp(os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta.temp"))
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "darkmatter_{tool}_{sample}_orfs.log")
  conda:
    PALM
  shell:
    """
    # Ativa modo de debug para ver o erro no log
    set -x 
    
    # Limpa arquivo de saida
    > {output.orfs}

    # Evita falha do getorf caso o arquivo de entrada nao possua sequencias
    if [ ! -s {input.fasta} ]; then
        echo "[INFO] Nenhum nohit identificado. Gerando arquivo vazio." >> {log}
        exit 0
    fi

    # Lista explicita para evitar erro de expansao do shell
    for table in 1 3 4 5 6 11 16; do
        tmp_file="$(dirname {output.orfs})/tmp_Orfs_{wildcards.sample}_{wildcards.tool}_Table${{table}}.fasta"     
        getorf -sequence {input.fasta} -minsize 600 -table "$table" -find 0 -outseq "$tmp_file" 2>> {log} || true
        if [ -s "$tmp_file" ]; then
            sed -i.bak "s/^>/>gc_${{table}}_/g" "$tmp_file" 2>> {log}
            cat "$tmp_file" >> {output.orfs}
        fi
        
        rm -f "$tmp_file" "$tmp_file.bak"
    done
    """

#######################################################
# --- 3. Remove Redundancy from final Orfs Fasta --- #
###################################################### 
rule cd_hit:
  """
  Remove duplicated ORFs using CD-HIT.
  """
  input:
    fasta = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta.temp")
  output:  
    orfs = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "darkmatter_{tool}_{sample}_cdhit.log")
  conda:
    PALM
  shell:
    """
    if [ -s {input.fasta} ]; then
        #remove duplicated orfs
        cd-hit -i {input.fasta} -o {output.orfs} -c 0.9 -d 1 -T {resources.threads} -M {resources.mem_mb} 2>> {log} 2>&1

        #substitute spaces per _
        sed -i.bak '/^>/ s/ /_/g' {output.orfs} >> {log} 2>&1
        rm -f {output.orfs}.bak
    else
        echo "[INFO] Nenhuma sequência para agrupar. Criando arquivo vazio." > {log}
        touch {output.orfs}
    fi
    """

###################################################
# --- 4. Iddentify putative RdRp within Data --- #
################################################## 
rule palm_annot:
  input:
    fasta = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta")
  output:
    fev  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp.fev"),
    rdrp = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp.fasta")
  params:
    seqtype = config["palm_annot"]["seqtype"],
    minscore = config["palm_annot"]["minscore"],
    minpssmscore = config["palm_annot"]["minpssmscore"],
    palm_annot_dir = config["resources"]["palm_annot_dir"]
  conda:
    PALM
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "darkmatter_{tool}_{sample}_palmannot.log")  
  script:
    "../scripts/palm_annot_run.py"
  
##################################
# --- 5. Convert FEV to TSV --- #
################################# 

rule fev2tsv_single:
  input:
    fev = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp.fev")
  output:
    tsv = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp.tsv")
  params:
    palm_annot_dir = config["resources"]["palm_annot_dir"]
  conda:
    PALM
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "darkmatter_{tool}_{sample}_fev2tsv.log")
  script:
    "../scripts/fev2tsv_run.py"

##############################################
# --- 6. Summarize Dusk Matter Findings --- #
############################################## 
rule report_summarize:
    input:
        fasta = get_contigs_path,
        lineages = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_lineage.tsv"),
        #basta_lca = os.path.join(OUT_DIR, "{sample}", "basta_{tool}", "{sample}_{tool}_lca.tsv"),
        dusk = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp.tsv")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_Report_Diversity.html")
    params:
        out_dir = directory(os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}")),
        logos = "resources/logo/"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_report_summarize_{tool}.log")
    conda:
        REPORT
    shell:
        """
        # Passamos o arquivo do BASTA no argumento --diamond_path.
        # O script R foi atualizado para detectar automaticamente que é um arquivo BASTA.
        Rscript workflow/scripts/generate_report.R \
            --sample_name {wildcards.sample} \
            --fasta_path {input.fasta} \
            --diamond_path {input.lineages} \
            --darkmatter_path {input.dusk} \
            --output_dir {params.out_dir} \
            --report_name $(basename {output.html}) \
            --logos {params.logos} \
            --input workflow/scripts/Report_Model.Rmd \
            > {log} 2>&1
        """

#########################################################
# --- 7. Extract ORFs Putative Raw RdRp and Contigs --- #
#########################################################
rule dm_validate:
  input:
    fasta = get_contigs_path,
    orfs  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta"),
    dusk  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp.tsv")
  output:
    raw_rdrp    = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp_Orfs.fasta"),
    raw_contigs = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_RdRp_contigs.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "{sample}_dm_validate_{tool}.log")
  conda:
    REPORT
  run:
    import pandas as pd
    import re
    from Bio import SeqIO

    with open(log[0], "w") as logf:
        logf.write(f"[INFO] dm_validate started for {wildcards.sample} / {wildcards.tool}\n")

        # Lê o TSV — pode estar vazio se nenhum RdRp foi encontrado
        try:
            df = pd.read_csv(input.dusk, sep='\t')
        except Exception as e:
            logf.write(f"[ERROR] Could not read RdRp TSV: {e}\n")
            df = pd.DataFrame()

        if df.empty or 'Label' not in df.columns:
            logf.write("[WARNING] No RdRp hits found. Creating empty output files.\n")
            open(output.raw_rdrp, 'w').close()
            open(output.raw_contigs, 'w').close()
        else:
            target_orfs = set(df['Label'].astype(str).tolist())
            target_contigs = set()

            # Lógica para extrair o ID do Contig a partir do Label do ORF
            # Exemplo Label: gc_1_NODE_244_length_2206_cov_6.211439_1_[2204_-_3]_(REVERSE_SENSE)_
            # Queremos extrair: NODE_244_length_2206_cov_6.211439
            for label in target_orfs:
                match = re.search(r'^gc_\d+_(.+?)_\d+_\[', label)
                if match:
                    contig_id = match.group(1)
                    target_contigs.add(contig_id)
                else:
                    logf.write(f"[WARNING] Não foi possível parsear o contig ID do label: {label}\n")

            # Salva os ORFs identificados
            with open(output.raw_rdrp, "w") as out_orf:
                for record in SeqIO.parse(input.orfs, "fasta"):
                    if record.id in target_orfs:
                        SeqIO.write(record, out_orf, "fasta")

            # Salva os Contigs originais correspondentes
            with open(output.raw_contigs, "w") as out_contig:
                for record in SeqIO.parse(input.fasta, "fasta"):
                    if record.id in target_contigs:
                        SeqIO.write(record, out_contig, "fasta")

            logf.write(f"[INFO] Done. {len(target_orfs)} RdRp ORFs, {len(target_contigs)} contigs extracted.\n")