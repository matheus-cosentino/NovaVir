#!/usr/bin/env python3
"""
Generate Kraken2 classification statistics for MultiQC
Extracts classification counts and summary information
"""

import json
import sys

def main():
    try:
        kraken_report = snakemake.input.report
        output_json = snakemake.output.json
        output_tsv = snakemake.output.tsv
        sample_name = snakemake.wildcards.sample
        tool_name = snakemake.wildcards.tool
        log_file = snakemake.log[0]
    except NameError:
        print("ERROR: This script is designed to be run via Snakemake.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Could not access Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    log = open(log_file, 'w')
    
    try:
        # Parse Kraken2 report format:
        # %    clade_count    direct_count    rank_code    tax_id    scientific_name
        
        total_reads = 0
        classified_reads = 0
        unclassified_reads = 0
        classified_pct = 0.0
        
        with open(kraken_report, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 6:
                    pct = float(parts[0])
                    clade_count = int(parts[1])
                    rank_code = parts[3]
                    sci_name = parts[5].strip()
                    
                    # Get totals from the root (U = unclassified, no rank_code)
                    if rank_code == 'U':
                        unclassified_reads = clade_count
                    elif rank_code == '-' and 'Bacteria' not in sci_name and 'Archaea' not in sci_name and 'Eukaryota' not in sci_name:
                        # Root node (no rank) gives us total
                        total_reads = clade_count
                    elif rank_code == '-':
                        # Top-level domain
                        if total_reads == 0:
                            total_reads = clade_count
        
        # Calculate classified
        if total_reads > 0:
            classified_reads = total_reads - unclassified_reads
            classified_pct = (classified_reads / total_reads) * 100
        
        # Create statistics dictionary
        stats_dict = {
            "sample": sample_name,
            "tool": tool_name,
            "total_reads": int(total_reads),
            "classified_reads": int(classified_reads),
            "unclassified_reads": int(unclassified_reads),
            "classified_percent": round(classified_pct, 2)
        }
        
        # Write JSON
        with open(output_json, 'w') as f:
            json.dump(stats_dict, f, indent=2)
        
        # Write TSV
        with open(output_tsv, 'w') as f:
            f.write("Metric\tValue\n")
            f.write(f"Total Sequences\t{total_reads}\n")
            f.write(f"Classified Sequences\t{classified_reads}\n")
            f.write(f"Unclassified Sequences\t{unclassified_reads}\n")
            f.write(f"Classification Rate (%)\t{classified_pct:.2f}\n")
        
        log.write(f"Successfully generated Kraken2 statistics for {sample_name} ({tool_name})\n")
        log.write(f"Total sequences: {total_reads}, Classified: {classified_reads}, Classification rate: {classified_pct:.2f}%\n")
        
    except Exception as e:
        log.write(f"ERROR: {str(e)}\n")
        import traceback
        log.write(traceback.format_exc())
        log.close()
        sys.exit(1)
    
    log.close()

if __name__ == "__main__":
    main()
