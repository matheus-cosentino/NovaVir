#MSc. Matheus Cosentino


## imports
import os, re, glob, time, sys, subprocess, platform, yaml
from snakemake.io import expand
from collections import defaultdict

#collors
# collor palletes to be used
blue="\033[1;34m"  # blue
green="\033[1;32m" # green
red="\033[1;31m"   # red
ylo="\033[1;33m"   # yellow
nc="\033[0m"       # no color


## Functions

def get_final_outputs():
  final_outputs = []
  if MODULES["discovery"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_hits_with_lineage.tsv",
     sample = SAMPLE, source = SOURCE))
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/report/{{sample}}_Report_Diversity.html",
     sample = SAMPLE))
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.fasta",
     sample = SAMPLE))
  if MODULES["exploratory"]:  
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_{{source}}_hits.tsv",
     sample = SAMPLE, source = SOURCE))
  if MODULES["overview"]:
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/kraken2/{{sample}}_{{source}}_hits.tsv",
     sample = SAMPLE, source = SOURCE))
    final_outputs.append(expand(f"{config['output_dir']}/{{sample}}/trimmed/{{sample}}.html",
     sample = SAMPLE))

  return final_outputs

def define_input(data_dir):
    """
    Identifies and categorizes input files: contigs (.fasta),
    paired FastQ (*_R1.fastq.gz and *_R2.fastq.gz), and
    unpaired FastQ (*.fastq.gz without a recognized pair).

    Args:
        data_dir (str): The directory to search for input files.

    Returns:
        tuple: A tuple containing three dictionaries:
               (contigs, paired_fastq, unpaired_fastq), and a list of warnings.
    """
    all_files = glob.glob(os.path.join(data_dir, "*"))
    
    # Dictionaries to store categorized files
    contigs = {} # {sample_name: 'path/to/file.fasta'}
    fastq_mates = defaultdict(dict) # Temporary storage for FastQ pairs {sample: {1: R1_path, 2: R2_path}}
    unpaired_fastq = {} # {sample_name: 'path/to/file.fastq.gz'}
    
    warnings = []

    for file_path in all_files:
        basename = os.path.basename(file_path)
        
        # 1. Identify Contigs (.fasta)
        if basename.endswith(".fasta") or basename.endswith(".fa"):
            # Use the full filename as the sample ID for contigs
            sample = basename.rsplit('.', 1)[0]
            contigs[sample] = file_path
            
        # 2. Identify FastQ files (*.fastq.gz)
        elif basename.endswith(".fastq.gz"):
            # Try to match the paired-end naming pattern: SAMPLE_R[12].fastq.gz
            match = re.search(r"^(?P<sample>.+)_R(?P<mate>[12])\.fastq\.gz$", basename)
            
            if match:
                sample = match.group("sample")
                mate = match.group("mate")
                fastq_mates[sample][mate] = file_path
            else:
                # If it doesn't match the paired pattern, treat it as a potential unpaired file
                # Use the full filename minus the extension as the sample ID
                sample = basename.replace(".fastq.gz", "")
                unpaired_fastq[sample] = file_path
        
        # 3. Handle other files (optional)
        # else:
        #     warnings.append(f"{blue}[INFO]{nc} Ignoring file '{basename}'.")


    # 3. Process Paired/Unpaired FastQ
    paired_fastq = {}
    
    # Check the temporarily stored FastQ files for complete pairs
    for sample, mates in fastq_mates.items():
        if "1" in mates and "2" in mates:
            # Complete pair found: move R1 and R2 paths to the paired_fastq dictionary
            paired_fastq[sample] = mates
        else:
            # Incomplete pair: the single file found is treated as unpaired
            if "1" in mates:
                warnings.append(f"{red}[WARNING]{nc} Sample '{sample}' is incomplete (missing R2). R1 treated as unpaired.")
                unpaired_fastq[f"{sample}_R1"] = mates["1"]
            elif "2" in mates:
                warnings.append(f"{red}[WARNING]{nc} Sample '{sample}' is incomplete (missing R1). R2 treated as unpaired.")
                unpaired_fastq[f"{sample}_R2"] = mates["2"]
    
    
    # Summary Info
    warnings.append(f"\n{blue}[INFO]{nc} Total valid **Contig** samples: {len(contigs)}")
    warnings.append(f"{blue}[INFO]{nc} Total valid **Paired FastQ** samples: {len(paired_fastq)}")
    warnings.append(f"{blue}[INFO]{nc} Total valid **Unpaired FastQ** samples: {len(unpaired_fastq)}")

    return contigs, paired_fastq, unpaired_fastq, warnings
