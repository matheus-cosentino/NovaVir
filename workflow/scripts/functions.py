###################################################################################
# Snakefile
#MSc. Matheus Cosentino 
###################################################################################
# oooooooooo.    o8o                               oooooo     oooo  o8o           #
# `888'   `Y8b   `"'                                `888.     .8'   `"'           #
#  888      888 oooo   .oooo.o  .ooooo.   .ooooo.    `888.   .8'   oooo  oooo d8b #
#  888      888 `888  d88(  "8 d88' `"Y8 d88' `88b    `888. .8'    `888  `888""8P #
#  888      888  888  `"Y88b.  888       888   888     `888.8'      888   888     #
#  888     d88'  888  o.  )88b 888   .o8 888   888      `888'       888   888     #
# o888bood8P'   o888o 8""888P' `Y8bod8P' `Y8bod8P'       `8'       o888o d888b    #

###################################################################################
# version: 12.2025
###################################################################################

# collor palletes
blue="\033[1;34m"  # blue
green="\033[1;32m" # green
red="\033[1;31m"   # red
ylow="\033[1;33m"   # yellow
nc="\033[0m"       # no color                                                                               

###########################################
# --- 1. Import libraryes to be used --- #
##########################################
import os, re, glob, time, sys, subprocess, platform, yaml
from snakemake.io import expand
from collections import defaultdict
from Bio import SeqIO
import sys
import gzip

##########################################################
# --- 2. Function to obtain final outputs per module --- #
##########################################################
## in process
def get_final_outputs():
  final_outputs = []
  if MODULES["quality"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.json",
     sample = SAMPLE))


  if MODULES["reads"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.json",
     sample = SAMPLE))

  if MODULES["assembly"]:  
    final_outputs.append(,
     sample = SAMPLE, mate = MATE)

  if MODULES["kraken2"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/kraken2/{{sample}}_{{source}}_hits.tsv",
     sample = SAMPLE, source = SOURCE))

  if MODULES["diamond"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_hits_with_lineage.tsv",
     sample = SAMPLE, source = SOURCE))

  if MODULES["duskmatter"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/kraken2/{{sample}}_{{source}}_hits.tsv",
     sample = SAMPLE, source = "contigs"))
    
  return final_outputs



#############################################################
# --- 3. Clean Up SRA Downloads after total processing --- #
############################################################
## ready to test
def cleanup_downloaded_fastqs(log):
    """
    Delete Fastq Data downloaded from SRA after workflown conclusion,
    unless 'keep_download' is TRUE.
    """
    options = config.get("modules", {})
    download_only = options.get("keep_download", False)

    # Verify if there is any SRA within metadata
    sra_samples_exist = any(meta['mode'] == 'SRA' for meta in SAMPLE_META.values())

    # 1: Not possible to clean (No SRA)
    if not sra_samples_exist:
        return

    # 2: User delanded to keep download
    if download_only:
        print(f"\n{blue}[INFO]{nc} Workflow (download_only) complete. *Keeping* SRA downloaded FASTQ files.")
        return

    # 3: Cleaning Space After Workflow completion
    print(f"\n{blue}[INFO]{nc} --- INITIATING POST-SUCCESS CLEANUP ---")
    print("Workflow complete. Removing temporary SRA FASTQ files...")
    
    files_removed = 0
    
    # Iteration over mmeta global in SAMPLE_META
    for sample, meta in SAMPLE_META.items():
        # Clean only SRA flag
        if meta['mode'] == 'SRA':
            # The list possesses paths
            for file_path in meta['files']:
                if os.path.exists(file_path):
                    try:
                        os.remove(file_path)
                        print(f"Removed: {file_path}")
                        files_removed += 1
                    except OSError as e:
                        print(f"{red}[WARNING]{nc} Could not remove {file_path}: {e}")
            
    print(f"{blue}[INFO]{nc} Cleanup complete. {files_removed} files removed.")
    print("---------------------------------------\n")


###########################################################
# --- 4. Obtain Novel Fasta with no Hits Iddentified --- #
##########################################################
## ready to test
def filter_nohits():
    # Access snakemake variables
    try:
        blast_file = snakemake.input.diamond
        fasta_file = snakemake.input.fasta
        output_file = snakemake.output.nohits
        sample_name = snakemake.wildcards.sample
        log_file = snakemake.log[0]
    except NameError:
        print("ERROR: This script is designed to be run via Snakemake.", file=sys.stderr)
        print("It expects 'snakemake.input', 'snakemake.output', etc.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Could not access Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    # Configure logging to stderr (which snakemake redirects to the log file)
    print(f"Starting sample {sample_name}: filtering no-hit sequences", file=sys.stderr)
    print(f"  Log file: {log_file}", file=sys.stderr)
    print(f"  Input FASTA: {fasta_file}", file=sys.stderr)
    print(f"  DIAMOND results: {blast_file}", file=sys.stderr)
    print(f"  Output FASTA: {output_file}", file=sys.stderr)
 
    # Read BLAST hits
    ids_with_hits = set()
    try:
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.strip():
                    parts = line.split('\t')
                    if parts:  # Ensure line has content
                        ids_with_hits.add(parts[0])
    except FileNotFoundError:
        print(f"ERROR: DIAMOND hits file not found at {blast_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading DIAMOND file {blast_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(ids_with_hits)} unique sequences with hits", file=sys.stderr)
    
    # Process and write sequences in one pass
    count_no_hits = 0
    count_total = 0
    
    try:
        # Handle gzipped input if needed
        open_func = open
        read_mode = 'rt'
        if fasta_file.endswith('.gz'):
            open_func = gzip.open
            print("Detected gzipped FASTA input", file=sys.stderr)
            
        with open_func(fasta_file, read_mode) as fasta_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                count_total += 1
                if record.id not in ids_with_hits:
                    count_no_hits += 1
                    SeqIO.write(record, output_handle, "fasta")
                    
    except FileNotFoundError:
        print(f"ERROR: FASTA file not found at {fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR processing FASTA file {fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Summary statistics
    print(f"Processed {count_total} total sequences", file=sys.stderr)
    print(f"Saved {count_no_hits} sequences with no hits ({count_no_hits/max(1,count_total)*100:.1f}%)", file=sys.stderr)
    print(f"Output written to: {output_file}", file=sys.stderr)
    print("Filtering complete.", file=sys.stderr)


####################################################
# --- 5. Obtain Metadata to Library Processing --- #
####################################################
## ready to test
def identify_data_type(sample_list, data_dir):
    """
    Verify if each sample within the file is a existing contig file, paired fastq or unpaired fastq.
    If none is found, signals for download (SRA).
    
    Returns:
        sample_meta (dict): Dictionary with metadada for sample.
                            Structure: {
                                'sample_A': {'mode': 'PAIRED', 'files': [path1, path2]},
                                'sample_B': {'mode': 'SRA', 'files': [path1_future, path2_future]}
                            }
    """
    sample_meta = {}
    
    print(f"[INFO] Checking input availability for {len(sample_list)} samples...")

    for sample in sample_list:
        # 1. Definição de caminhos esperados (Prioridade de checagem)
        
        # A. Contigs (.fasta or .fa or .fas)
        path_fasta = os.path.join(data_dir, f"{sample}.fasta")
        path_fa = os.path.join(data_dir, f"{sample}.fa")
        path_fas = os.path.join(data_dir, f"{sample}.fas")
        
        # B. Paired FastQ (_R1/_R2 ou _1/_2)
        p_r1 = os.path.join(data_dir, f"{sample}_R1.fastq.gz")
        p_r2 = os.path.join(data_dir, f"{sample}_R2.fastq.gz")
        p_1  = os.path.join(data_dir, f"{sample}_1.fastq.gz")
        p_2  = os.path.join(data_dir, f"{sample}_2.fastq.gz")
        
        # C. Unpaired FastQ
        p_unpaired = os.path.join(data_dir, f"{sample}.fastq.gz")

        # 1. Contig?
        if os.path.exists(path_fasta):
            sample_meta[sample] = {'mode': 'CONTIGS', 'files': [path_fasta]}
        elif os.path.exists(path_fa):
            sample_meta[sample] = {'mode': 'CONTIGS', 'files': [path_fa]}
        elif os.path.exists(path_fas):
            sample_meta[sample] = {'mode': 'CONTIGS', 'files': [path_fas]}
            
        # 2. Paired?
        elif os.path.exists(p_r1) and os.path.exists(p_r2):
            sample_meta[sample] = {'mode': 'PAIRED', 'files': [p_r1, p_r2]}
        elif os.path.exists(p_1) and os.path.exists(p_2):
            sample_meta[sample] = {'mode': 'PAIRED', 'files': [p_1, p_2]}
            
        # 3. Unpaired?
        elif os.path.exists(p_unpaired):
            sample_meta[sample] = {'mode': 'UNPAIRED', 'files': [p_unpaired]}
            
        # 4. None? Download (SRA)
        else:
            # Warning; Define the place to find file after the download for DAG Build
            future_r1 = os.path.join(data_dir, f"{sample}_1.fastq.gz")
            future_r2 = os.path.join(data_dir, f"{sample}_2.fastq.gz")
            
            sample_meta[sample] = {'mode': 'SRA', 'files': [future_r1, future_r2]}
            print(f"[WARNING] Input missing for '{sample}'. Marked for SRA Download.")

    return sample_meta


####################################################################
# --- 6. Obtain path to Library Reads availablle or downloaded --- #
####################################################################
## ready to test
def get_input_r1(wildcards):
    """Returns R1 local or future"""
    meta = SAMPLE_META[wildcards.sample]
    if meta['mode'] in ['PAIRED', 'SRA']:
        return meta['files'][0]
    if meta['mode'] == 'UNPAIRED':
        return meta['files'][0]
    return []

def get_input_r2(wildcards):
    """Returns R2 local or future"""
    meta = SAMPLE_META[wildcards.sample]
    if meta['mode'] in ['PAIRED', 'SRA']:
        return meta['files'][1]
    return []