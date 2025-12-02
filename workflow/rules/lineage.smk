# workflow/rules/lineage.smk

rule map_accession_to_taxid:
    """
    Maps protein IDs (Subject ID, column 2 of DIAMOND) to TaxIDs.
    """
    input:
        hit_file="{out_dir}/{sample}/{tool}/diamond/{sample}_contigs_report.txt",
        taxid_map="resources/database/prot.accession2taxid.gz"
    output:
        temp(f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_hits_with_taxid.tmp")
    shadow: "minimal" 
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_{{source}}_map_taxid.log"
    shell:
        """
        # Extract unique protein IDs (skip header if present)
        tail -n +2 {input.hit_file} | cut -f 2 | sort -u > {output}.protein_ids.tmp
        
        # Filter taxid map and create lookup
        zcat {input.taxid_map} | grep -Fwf {output}.protein_ids.tmp > {output}.filtered_map.tmp
        
        # Add taxid to hits
        awk -F'\\t' -v OFS='\\t' '
        NR==FNR {{
            # Store both accession and accession.version as keys
            taxid_map[$1] = $3
            taxid_map[$2] = $3
            next
        }}
        FNR==1 {{
            # Print header with added taxid column
            print $0, "taxid"
            next
        }}
        {{
            protein_id = $2
            taxid = (protein_id in taxid_map) ? taxid_map[protein_id] : "NOT_FOUND"
            print $0, taxid
        }}' {output}.filtered_map.tmp {input.hit_file} > {output} 2>> {log}
        
        # Cleanup
        rm {output}.protein_ids.tmp {output}.filtered_map.tmp
        """

rule split_hits_by_taxid:
    """
    Filters the input file to keep only hits with valid TaxIDs.
    """
    input:
        f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_hits_with_taxid.tmp"
    output:
        valid_hits=temp(f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_valid_hits.temp")
    params:
        header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    shell:
        """
        # 1. Create the output file with the header
        echo -e "{params.header}" > {output.valid_hits}
        
        # 2. Filter data: Skip header (NR==1) AND only print if last column is NOT "NOT_FOUND"
        awk -F'\\t' '
        NR==1 {{ next }} 
        $NF != "NOT_FOUND" {{
            print $0 >> "{output.valid_hits}"
        }}' {input}
        """


rule append_lineage:
    input:
        valid_hits=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_valid_hits.temp",
        nodes="resources/database/nodes.dmp",
        names="resources/database/names.dmp"
    output:
        f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_hits_with_lineage.tsv"
    params:
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tTaxid",
        lineage_header="Lineage\tCelular\tAcelular\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_{{source}}_get_lineage.log"
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {input.nodes})
        
        # Skip header and extract taxids
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # 1. Get lineage information, outputting TaxID (default col 1) + 12 reformatted columns (col 2-13), all separated by TAB
        # Usamos '\t' para garantir que os campos sejam separados por tabulação.
        taxonkit lineage --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \\
        taxonkit reformat --data-dir "${{DB_DIR}}" \\
            -f "{{C}}\\t{{a}}\\t{{d}}\\t{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" \\
             2>> {log} > {output}.lineage.tmp
        
        # 2. Extract ONLY the 12 reformatted columns (cols 2-13), dropping the TaxID (col 1),
        # as the TaxID is already available in the main DIAMOND output (col 13).
        # We use cut on the tab-separated intermediate file.
        cut -f 2- {output}.lineage.tmp > {output}.reformat_only.tmp
        
        # 3. Create final output with comprehensive header
        echo -e "{params.base_header}\\t{params.lineage_header}" > {output}
        
        # 4. Combine original DIAMOND data (cols 1-12) + Taxid (col 13) + Reformatted Ranks (cols 14-25)
        tail -n +2 {input.valid_hits} | cut -f 1-12 | \\
        paste - <(tail -n +2 {input.valid_hits} | cut -f 13) {output}.reformat_only.tmp >> {output}
        
        # Cleanup
        rm {output}.taxids.tmp {output}.lineage.tmp {output}.reformat_only.tmp
        """
