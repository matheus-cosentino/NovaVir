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

rule rvdb_validate_structural:
  input:
    fasta = get_contigs_path,
    orfs  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta"),
    tsv   = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_RVDB_Summary.csv")
  output:
    raw_rvdb_orfs    = os.path.join(OUT_DIR, "{sample}", "rvdb_structural_{tool}", "{sample}_structural_Orfs.fasta"),
    raw_rvdb_contigs = os.path.join(OUT_DIR, "{sample}", "rvdb_structural_{tool}", "{sample}_structural_contigs.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "{sample}_rvdb_validate_structural_{tool}.log")
  conda:
    REPORT
  run:
    import pandas as pd
    import re
    from Bio import SeqIO
    
    # Carrega o CSV sumarizado do RVDB
    try:
        df = pd.read_csv(input.tsv, sep='\t')
        
        # --- FILTRO DE CONFIANÇA ALTA APENAS ---
        if 'Confidence' in df.columns:
            df = df[df['Confidence'] == 'High']
            
        target_orfs = set(df['Sequence_ID'].astype(str).tolist())
    except (pd.errors.EmptyDataError, KeyError):
        target_orfs = set()
        
    target_contigs = set()

    # Logica para extrair o ID do Contig a partir do Label do ORF
    for label in target_orfs:
        match = re.search(r'^gc_\d+_(.+?)_\d+_\[', label)
        if match:
            contig_id = match.group(1)
            target_contigs.add(contig_id)
        else:
            with open(log[0], "a") as f:
                f.write(f"[WARNING] Nao foi possivel parsear o contig ID do label: {label}\n")

    # Salva os ORFs identificados
    with open(output.raw_rvdb_orfs, "w") as out_orf:
        if len(target_orfs) > 0:
            for record in SeqIO.parse(input.orfs, "fasta"):
                if record.id in target_orfs:
                    SeqIO.write(record, out_orf, "fasta")

    # Salva os Contigs originais correspondentes 
    with open(output.raw_rvdb_contigs, "w") as out_contig:
        if len(target_contigs) > 0:
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id in target_contigs:
                    SeqIO.write(record, out_contig, "fasta")

rule rvdb_validate_pol:
  input:
    fasta = get_contigs_path,
    orfs  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_ORFs.fasta"),
    tsv   = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_RVDB_Summary.csv")
  output:
    raw_rvdb_orfs    = os.path.join(OUT_DIR, "{sample}", "rvdb_pol_{tool}", "{sample}_pol_Orfs.fasta"),
    raw_rvdb_contigs = os.path.join(OUT_DIR, "{sample}", "rvdb_pol_{tool}", "{sample}_pol_contigs.fasta")
  log:
    os.path.join(OUT_DIR, "{sample}", "log", "{sample}_rvdb_validate_pol_{tool}.log")
  conda:
    REPORT
  run:
    import pandas as pd
    import re
    from Bio import SeqIO
    
    # 1. Carrega o CSV tratando as aspas que aparecem nos seus dados 
    try:
        # O quotechar='"' garante que o pandas remova as aspas dos IDs e anotações 
        df = pd.read_csv(input.tsv, sep='\t', quotechar='"')
        
        # 2. Filtro atualizado para incluir "reverse" (comum no RVDB para pols) 
        pol_pattern = r'pol|polymerase|transcriptase|polyprotein|rdrp|integrase|protease|retrotransposon|reverse|\brt\b'
        
        mask = pd.Series(False, index=df.index)
        for col in ['Annotation', 'Description', 'str']:
            if col in df.columns:
                mask |= df[col].str.contains(pol_pattern, case=False, na=False)
        
        df = df[mask]
        target_orfs = set(df['Sequence_ID'].astype(str).tolist())
    except (pd.errors.EmptyDataError, KeyError):
        target_orfs = set()
        
    target_contigs = set()

    # 3. Extração do ID do Contig via Regex 
    for label in target_orfs:
        # Captura o que está entre 'gc_X_' e o penúltimo '_X_[' 
        match = re.search(r'gc_\d+_(.+?)_\d+_\[', label)
        if match:
            contig_id = match.group(1)
            target_contigs.add(contig_id)
        else:
            with open(log[0], "a") as f:
                f.write(f"[WARNING] Falha no parse do ID: {label}\n") 

    # 4. Salva os ORFs 
    with open(output.raw_rvdb_orfs, "w") as out_orf:
        if len(target_orfs) > 0:
            for record in SeqIO.parse(input.orfs, "fasta"):
                if record.id in target_orfs: 
                    SeqIO.write(record, out_orf, "fasta") 

    # 5. Salva os Contigs (usando record.id completo) 
    with open(output.raw_rvdb_contigs, "w") as out_contig:
        if len(target_contigs) > 0:
            for record in SeqIO.parse(input.fasta, "fasta"):
                # Se o seu FASTA tem o nome completo (NODE_10426_length...), 
                # o record.id deve bater com o contig_id extraído.
                if record.id in target_contigs: 
                    SeqIO.write(record, out_contig, "fasta") 