# workflow/rules/duskmatter.smk

rule get_nohit_fasta:
  """
    Get contigs within no hits against the nr Diamond database
  """
  input:
    fasta=f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/contigs.fasta",
    diamond=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_with_lineage.tsv"
  output:
    nohits=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_contigs_nohit.fasta"
  log:
    f"{config['output_dir']}/{{sample}}/{{sample}}_nohits.log"
  conda:
    SCRIPTS
  script:
    "scripts/filter_nohits.py"


rule find_orfs:
  """
    Extract orfs in distinct genetic codes used by RNA Viruses.
  """
  input: 
    fasta=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_contigs_nohit.fasta"
  output:
    orfs=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_ORFs.fasta.tmp"
  log:
    f"{config['output_dir']}/{{sample}}/logs/{{sample}}_orfs.log"
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

rule cd_hit:
  """
  Remove duplicated ORFs using CD-HIT.
  """
  input:
    fasta=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_ORFs.fasta.tmp"
  output:
    orfs=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_ORFs.fasta"
  log:
    f"{config['output_dir']}/{{sample}}/logs/{{sample}}_cdhit.log"
  conda:
    PALM
  shell:
    """
    #remove duplicated orfs
    cd-hit -i {input.fasta} -o {output.orfs} -c 0.98 -d 1 -T {resources.threads} >> {log} 2>&1

    #substitute spaces per _
    sed -i.bak '/^>/ s/ /_/g' {output.orfs} >> {log} 2>&1
    rm -f {output.orfs}.bak
    """

rule palm_annot:
  input:
    fasta = f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_ORFs.fasta"
  output:
    fev=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.fev",
    rdrp = f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.fasta"
  params:
    seqtype = config["params"]["palm_annot"]["seqtype"],
    minscore = config["params"]["palm_annot"]["minscore"],
    minpssmscore = config["params"]["palm_annot"]["minpssmscore"],
    palm_annot_dir = config["db"]["palm_annot_dir"]
    palm_annot_script = f"{config['db']['palm_annot_dir']}/py/palm_annot.py"
  conda:
    PALM
  log:
    f"{config['output_dir']}/{{sample}}/logs/{{sample}}_palmannot.log"
  shell:
    """
    export PATH={params.palm_annot_dir}/bin:{params.palm_annot_dir}/py:$$PATH
        
    python3 {params.palm_annot_script} \
     --input {input.fasta} \
     --seqtype {params.seqtype} \
     --fev {output.fev} \
     --rdrp {output.rdrp} \
     --minscore {params.minscore} \
     --threads {resources.threads} \
     --minpssmscore {params.minpssmscore} \
      2> {log}
        """

# Regra para converter FEV para TSV 
rule fev2tsv_single:
  input:
    fev=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.fev"
  output:
    tsv=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.tsv"
  params:
    palm_annot_dir = config["db"]["palm_annot_dir"]
  conda:
    PALM
  log:
    f"{config['output_dir']}/{{sample}}/logs/fev2tsv_{{sample}}.log"
  shell:
    """
    export PATH={params.palm_annot_dir}/bin:{params.palm_annot_dir}/py:$$PATH
    fev2tsv.py --input {input.fev} --output {output.tsv} 2> {log}
    """

