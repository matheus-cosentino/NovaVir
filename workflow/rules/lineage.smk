###################################################################################
#                       workflow/rules/lineage.smk                                #
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


#################################################
# --- 1. Map Proteins Id Diamond for Taxid --- #
################################################ 

rule map_accession_to_taxid:
    """
    Maps protein IDs (Subject ID, column 2 of DIAMOND) to TaxIDs.
    """
    input:
        hit_file = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_report.txt")
    output:
        ids = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_taxid.tmp")
    shadow: 
        "minimal" 
    params:
        #taxid_map="resources/database/prot.accession2taxid.gz"
        taxid_map = config["resources"]["taxonmap"]
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_map_acc_prot.log")
    shell:    
       """
        
        echo "Starting taxid mapping with join approach..." > {log}
        
        # 1. Extract protein IDs
        HEADER=$(head -1 {input.hit_file} | grep -c "^qseqid" || echo "0") # [cite: 24]
        
        echo "Extract step of protein ID Header..." >> {log}
        
        # ... (Protein ID extraction logic remains unchanged, around lines 25-29) ...
        if [ "$HEADER" -eq 1 ];
        then
            tail -n +2 {input.hit_file} | cut -f 2 | sort -u > {output}.protein_ids.tmp
            echo "Header idd in the diamond file..." >> {log}
        else
            cut -f 2 {input.hit_file} | sort -u > {output}.protein_ids.tmp
            echo "Header idd NOT the diamond file..." >> {log} # [cite: 26]
        fi

        echo "step of protein ID finished ..." >> {log}
        
        # Check if we have protein IDs
        PROTEIN_COUNT=$(wc -l < {output}.protein_ids.tmp) # [cite: 27]
        if [ "$PROTEIN_COUNT" -eq 0 ];
        then
            echo "WARNING: No protein IDs found. Creating empty output with header." >> {log} # [cite: 28, 29]
            echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\ttaxid" > {output}
            exit 0
        fi
        
        # Add tab for joining on second column
        awk '{{print $1 "\\t"}}' {output}.protein_ids.tmp > {output}.protein_ids_for_join.tmp
        echo "Prepared {output}.protein_ids_for_join.tmp with $PROTEIN_COUNT unique protein IDs." >> {log}
        
        # 2. Check if taxid map is gzipped and define $DECOMPRESS
        DECOMPRESS="cat" # Default to cat 
        if [[ {params.taxid_map} == *.gz ]];
        then
            DECOMPRESS="zcat"
            echo "Taxid map file is *.gz, using zcat." >> {log}
        fi
        
        # 3. Check taxid map header and prepare it
        # Temporarily redirect errors from $DECOMPRESS to a side file for debugging
        echo "Attempting to read header with $DECOMPRESS command..." >> {log}
        
        # Capture error and check exit code for decompression failure
        $DECOMPRESS {params.taxid_map} 2> {log}.decompress_err | head -1 | grep -c "^accession" > {log}.header_check.tmp
        TAXID_HEADER=$(cat {log}.header_check.tmp || echo "0")
        
        if [ $? -ne 0 ]; then
             echo "FATAL ERROR: Decompression/read of taxid map failed." >> {log}
             cat {log}.decompress_err >> {log}
             exit 1
        fi
        
        echo "Taxid map header present: $TAXID_HEADER" >> {log}
        
        # 4. Process the taxid map file and sort it
        echo "Extracting columns from taxid map..." >> {log} # [cite: 33, 34, 35]
        
        if [ "$TAXID_HEADER" -eq 1 ];
        then
            # Skip header and extract columns 2 (accession.version) and 3 (taxid)
            echo "Extracting columns from taxid map (skipping header)..." >> {log}
            $DECOMPRESS {params.taxid_map} | tail -n +2 | cut -f 2,3 | sort -k1,1 > {output}.taxmap_sorted.tmp
        else
            # Extract columns 2 and 3 directly
            echo "Extracting columns from taxid map (no header)..." >> {log}
            $DECOMPRESS {params.taxid_map} | cut -f 2,3 | sort -k1,1 > {output}.taxmap_sorted.tmp
        fi
        
        # --- NEW DEBUG LINE ---
        TAXMAP_LINES=$(wc -l < {output}.taxmap_sorted.tmp)
        echo "Taxid map sorted with $TAXMAP_LINES lines" >> {log} # This must be written now!
        
        # 5. Join on first field (accession.version)
        echo "Joining protein IDs with taxid map..." >> {log} # [cite: 36]
        join -t $'\\t' -1 1 -2 1 -a 1 -e "NOT_FOUND" -o 1.1,2.2 \
            <(sort -k1,1 {output}.protein_ids_for_join.tmp) \
            {output}.taxmap_sorted.tmp > {output}.id_to_taxid.tmp 2>> {log} # [cite: 37]
        
        # ... (rest of the rule remains unchanged, around lines 37-44) ...

        JOIN_RESULT=$?
        JOIN_LINES=$(wc -l < {output}.id_to_taxid.tmp 2>/dev/null || echo "0")
        echo "Join completed with exit code $JOIN_RESULT, found $JOIN_LINES matches" >> {log} # [cite: 37]
        
        if [ $JOIN_RESULT -ne 0 ] || [ "$JOIN_LINES" -eq 0 ]; then
            echo "WARNING: Join may have failed or found no matches. Creating mapping with NOT_FOUND." >> {log} # [cite: 38, 39]
            # Create a mapping file with all NOT_FOUND
            awk '{{print $1 "\\tNOT_FOUND"}}' {output}.protein_ids.tmp > {output}.id_to_taxid.tmp
        fi
        
        # 6. Create a mapping file for awk
        awk 'BEGIN {{FS=OFS="\\t"}} {{print $1, $2}}' {output}.id_to_taxid.tmp > {output}.mapping.tmp # [cite: 40]
        
        # 7. Add taxid column to original file
        echo "Adding taxid column to original hits..." >> {log}
        awk 'BEGIN {{
            FS=OFS="\\t"
            while (getline < "{output}.mapping.tmp") {{
                taxid_map[$1] = $2
            }}
        }}
        {{
            if (NR == 1 && $1 == "qseqid") {{
                print $0, "taxid"
            }} else {{
                protein_id = $2
                taxid = (protein_id in taxid_map) ? taxid_map[protein_id] : "NOT_FOUND" # [cite: 41, 42]
                print $0, taxid
            }}
        }}' {input.hit_file} > {output} 2>> {log}
        
        AWK_RESULT=$?
        if [ $AWK_RESULT -ne 0 ]; then
            echo "ERROR: awk failed with exit code $AWK_RESULT" >> {log} # [cite: 43]
            exit 1
        fi
        
        # Cleanup
        rm -f {output}.*.tmp {log}.decompress_err {log}.header_check.tmp
        
        echo "Done." >> {log} # [cite: 44] 
       """


##################################################
# --- 2. Split Hits for Taxid and Not Found --- #
################################################# 

rule split_hits_by_taxid:
    """
    Filters the input file to keep only hits with valid TaxIDs.
    """
    input:
      os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_taxid.tmp")
    output:
      valid_hits=temp(os.path.join(OUT_DIR, "{sample}", "diamond_{source}",  "{sample}_{source}_valid_hits.tmp"))    
    params:
      header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    log:
      os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_split_hits.log")
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

################################################
# --- 3. Append taxonomic to Diamond file --- #
############################################### 

rule append_lineage:
    input:
       valid_hits = os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_valid_hits.tmp")    
    output:
       os.path.join(OUT_DIR, "{sample}", "diamond_{source}", "{sample}_{source}_hits_with_lineage.tsv")
    params:
        nodes = config["resources"]["taxonnodes"],
        names = config["resources"]["taxonnames"],
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tTaxid",
        lineage_header="Lineage\tCelular\tAcelular\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{source}_append_lineage.log")
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {params.nodes})
        
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
