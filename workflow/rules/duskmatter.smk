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
#                              version: 12.2025                                   #
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
    nohits = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_{tool}_nohit.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "duskmatter_{tool}_{sample}_nohits.log")
  conda:
    CORE
  script:
    "../scripts/filter_nohits.py"


#################################################
# --- 2. Find Orfs for the no Hits Contigs --- #
################################################ 
rule find_orfs:
  """
    Extract orfs in distinct genetic codes used by RNA Viruses.
  """
  input: 
    fasta = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_{tool}_nohit.fasta")
  output:
    orfs = temp(os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_ORFs.fasta.temp"))
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "duskmatter_{tool}_{sample}_orfs.log")
  conda:
    PALM
  shell:
    """
    genetic_codes="1 3 4 5 6 11 16"
    # Loop through each genetic code table and run getorf generating temp files

    for table in $genetic_codes; do
        getorf {input.fasta} -minsize 600 -table "$table" -find 0 -outseq "Orfs_{wildcards.sample}_Table${{table}}.fasta"
        sed -i "s/^>/>gc_${{table}}_/g" "Orfs_{wildcards.sample}_Table${{table}}.fasta"
    done
    
    # Combine orfs files
    cat Orfs_{wildcards.sample}_Table*.fasta > {output.orfs}  2>> {log}
    # Clean up temp files
    rm Orfs_{wildcards.sample}_Table*.fasta
    """

#######################################################
# --- 3. Remove Redundancy from final Orfs Fasta --- #
###################################################### 
rule cd_hit:
  """
  Remove duplicated ORFs using CD-HIT.
  """
  input:
    fasta = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_ORFs.fasta.temp")
  output:  
    orfs = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_ORFs.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "duskmatter_{tool}_{sample}_cdhit.log")
  conda:
    PALM
  shell:
    """
    #remove duplicated orfs
    cd-hit -i {input.fasta} -o {output.orfs} -c 0.9 -d 1 -T {resources.threads} >> {log} 2>&1

    #substitute spaces per _
    sed -i.bak '/^>/ s/ /_/g' {output.orfs} >> {log} 2>&1
    rm -f {output.orfs}.bak
    """

###################################################
# --- 4. Iddentify putative RdRp within Data --- #
################################################## 
rule palm_annot:
  input:
    fasta = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_ORFs.fasta")
  output:
    fev  = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_RdRp.fev"),
    rdrp = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_RdRp.fasta")
  params:
    seqtype = config["palm_annot"]["seqtype"],
    minscore = config["palm_annot"]["minscore"],
    minpssmscore = config["palm_annot"]["minpssmscore"],
    palm_annot_dir = config["resources"]["palm_annot_dir"]
  conda:
    PALM
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "duskmatter_{tool}_{sample}_palmannot.log")  
  script:
    "../scripts/palm_annot_run.py"
  
##################################
# --- 5. Convert FEV to TSV --- #
################################# 

rule fev2tsv_single:
  input:
    fev = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_RdRp.fev")
  output:
    tsv = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_RdRp.tsv")
  params:
    palm_annot_dir = config["resources"]["palm_annot_dir"]
  conda:
    PALM
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "duskmatter_{tool}_{sample}_fev2tsv.log")
  script:
    "../scripts/fev2tsv_run.py"

##############################################
# --- 6. Summarize Dusk Matter Findings --- #
############################################## 
rule report_summarize:
  """
  Executes an R script to generate a primary HTML report.
  """
  # ADICIONAR ESTA CONSTRAINT TAMBÉM
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    fasta = get_contigs_path,
    diamond = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_lineage.tsv"),
    duskmatter = os.path.join(OUT_DIR, "{sample}", "duskmatter_{tool}", "{sample}_RdRp.tsv")
  output:
    html = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_Report_Diversity.html"),
    contig_summary_tsv = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_{tool}_summary.tsv"),
    kingdom_tsv = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_{tool}_summary_per_Kingdom.tsv"),
    viral_diamond_tsv = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_Viral_Diamond.tsv"),
    viral_summary_tsv = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_viral_{tool}_summary.tsv"),
    kingdom_png = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_KINGDOM_Classification_Log10.png"),
    viral_png = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}", "{sample}_Virus_Classification_Log10.png")
  params:
    logo_dirs = config["resources"]["logo_dirs"],
    output_dir = os.path.join(OUT_DIR, "{sample}", "duskmatter_report_{tool}") + "/"
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "{sample}_report_summarize_{tool}.log")
  conda:
    REPORT
  shell:
    """
    # Execute the R script inside a subshell
      Rscript workflow/scripts/generate_report.R \
      --sample_name {wildcards.sample} \
      --fasta_path {input.fasta} \
      --diamond_path {input.diamond} \
      --duskmatter_path {input.duskmatter} \
      --output_dir {params.output_dir} \
      --report_name {wildcards.sample}_Report_Diversity.html \
      --logos {params.logo_dirs} \
      --input workflow/rules/scripts/Report_Model.Rmd \
      2>&1 > {log}
     """