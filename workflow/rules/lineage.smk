# workflow/rules/lineage.smk

rule map_accession_to_taxid:
    """
    Maps protein IDs (Subject ID, column 2 of DIAMOND) to TaxIDs.
    """
    input:
        hit_file=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits.tsv",
        taxid_map="resources/database/prot.accession2taxid.gz"
    output:
        f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_with_taxid.tmp"
    shadow: "minimal" 
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_contigs_map_taxid.log"
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
    Splits the input file into two: one with valid taxids and one with "NOT_FOUND".
    """
    input:
        f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_with_taxid.tmp"
    output:
        valid_hits=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_valid_hits.temp",
        no_lineage=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_no_lineage.temp"
    params:
        header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    shell:
        """
        # Create output files with headers
        echo -e "{params.header}" > {output.valid_hits}
        echo -e "{params.header}" > {output.no_lineage}
        
        # Split data (skip header from input)
        awk -F'\\t' '
        NR==1 {{ next }}  # Skip input header
        {{
            if ($NF == "NOT_FOUND") {{
                print $0 >> "{output.no_lineage}"
            }} else {{
                print $0 >> "{output.valid_hits}"
            }}
        }}' {input}
        """

rule append_lineage:
    """
    Takes the valid hits, gets the lineage for them, and creates the final output file.
    """
    input:
        valid_hits=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_valid_hits.temp",
        nodes="resources/database/nodes.dmp",
        names="resources/database/names.dmp"
    output:
        f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_with_lineage.tsv"
    params:
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore",
        lineage_header="Taxid\tLineage\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_contigs_get_lineage.log"
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {input.nodes})
        
        # Skip header and extract taxids
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # Get lineage information
        taxonkit lineage --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \
        taxonkit reformat --data-dir "${{DB_DIR}}" \
            -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F -t 2>> {log} > {output}.lineage.tmp
        
        # Create final output with comprehensive header
        echo -e "{params.base_header}\\t{params.lineage_header}" > {output}
        
        # Combine original data (without taxid header) with lineage data
        tail -n +2 {input.valid_hits} | cut -f 1-12 | \
        paste - <(tail -n +2 {input.valid_hits} | cut -f 13) <(cat {output}.lineage.tmp) >> {output}
        
        # Cleanup
        rm {output}.taxids.tmp {output}.lineage.tmp
        """