# Snakefile

import os
from glob import glob

# --- 1. Load Configuration ---
configfile: "config/workflow_config.yaml"

# --- 2. Discover Samples (MODIFIED FOR MIXED MODES) ---
RAW_DATA_DIR = "data/raw_data"
CONTIGS_DIR = "data/contigs"

# New: Dictionary to store the specific mode for each sample
SAMPLE_MODES = {}
# New: Flag to track if ANY sample requires SRA download (for cleanup logic)
ANY_SRA_DOWNLOADED = False 

# 1. Load Sample IDs from SRA file (which now serves as the master list)
try:
    sra_file_path = config["sra_file"]
    with open(sra_file_path, "r") as f:
        SAMPLES = [line.strip() for line in f if line.strip()]
    if not SAMPLES:
        raise ValueError("The sample list file is empty.")

except (KeyError, FileNotFoundError):
    raise ValueError(
        "CRITICAL: The 'sra_file' path is missing from config.yaml or the file is invalid. "
        "A list of sample IDs is required."
    )

# 2. Check the availability of each sample and assign its specific MODE
final_samples = []
for sample in SAMPLES:
    fastq_r1 = os.path.join(RAW_DATA_DIR, f"{sample}_1.fastq.gz")
    fastq_r2 = os.path.join(RAW_DATA_DIR, f"{sample}_2.fastq.gz")
    contig_glob_path = os.path.join(CONTIGS_DIR, f"{sample}*")

    # --- Mode Prioritization (Sample-specific assignment) ---

    # Priority 1: FASTQ Local (Full Pipeline) - Requires BOTH R1 and R2
    if os.path.exists(fastq_r1) and os.path.exists(fastq_r2):
        SAMPLE_MODES[sample] = 'FASTQ'
        
    # Priority 2: Contigs Pre-montados (Annotation Only Pipeline) - Check for .fasta or .fna
    elif glob(contig_glob_path + '.fasta') or glob(contig_glob_path + '.fna'):
        SAMPLE_MODES[sample] = 'CONTIGS'
        
    # Priority 3: SRA Download (Full Pipeline)
    else:
        SAMPLE_MODES[sample] = 'SRA'
        ANY_SRA_DOWNLOADED = True # Set flag if any sample needs SRA
        
    final_samples.append(sample)
        
SAMPLES = final_samples # Update SAMPLES list

print(f"Workflow Mode detected: MIXED (SRA, FASTQ, CONTIGS)")
print(f"Workflow will process the following samples:")
for sample, mode in SAMPLE_MODES.items():
    print(f"  {sample}: {mode}")

# --- 3. Define Global Environment Paths ---
DOWNLOAD = workflow.source_path(config["envs"]["download"])
FASTP = workflow.source_path(config["envs"]["fastp"])
SPADES = workflow.source_path(config["envs"]["spades"])
DIAMOND = workflow.source_path(config["envs"]["diamond"])
TAXONKIT = workflow.source_path(config["envs"]["taxonkit"])
SCRIPTS = workflow.source_path(config["envs"]["scripts"])
REPORT = workflow.source_path(config["envs"]["report"])
PALM = workflow.source_path(config["envs"]["palm_annot"])


# --- 4. Load the Workflow Rules ---
include: "workflow/rules/get_data.smk"
include: "workflow/rules/fastp.smk"
include: "workflow/rules/spades.smk"
include: "workflow/rules/diamond.smk"
include: "workflow/rules/lineage.smk"
include: "workflow/rules/reporting.smk"
include: "workflow/rules/duskmatter.smk"

# --- 4.5 Resolve Ambiguity (REQUIRED) ---
# Prioritize 'link_preassembled_contigs' over 'spades' for conditional fallback
#ruleorder: link_preassembled_contigs > spades


# --- 5. Obtain the final target (MODIFIED FOR MIXED MODES) ---
def get_final_targets(wildcards):
    final_files = []
    options = config.get("options", {})
    download_only = options.get("download_only", False)
    
    print(f"Workflow mode: FULL ANALYSIS / ANNOTATION (MIXED detected)")

    for sample in SAMPLES:
        mode = SAMPLE_MODES[sample]

        # 1. Modo de Download (Applies only if the sample's mode is SRA/FASTQ)
        if download_only and mode != 'CONTIGS':
            final_files.extend(
                expand(
                    [
                        "raw_data/{sample}_1.fastq.gz",
                        "raw_data/{sample}_2.fastq.gz"
                    ],
                    sample=[sample]
                )
            )
            continue # Move to the next sample

        # 2. Annotation Targets (Applies to all samples when not download_only)
        
        # Alvo principal de anotação de contigs
        final_files.append(
            f"{config['output_dir']}/{sample}/diamond/{sample}_contigs_hits_with_lineage.tsv"
        )
        
        # Alvo opcional de reads (ONLY for SRA/FASTQ modes)
        if options.get("diversity_reads", False) and mode != 'CONTIGS': 
            final_files.append(
                f"{config['output_dir']}/{sample}/diamond/{sample}_reads_hits_with_lineage.tsv"
            )
            
        # Alvo opcional de palm_annot (RdRp) (Applies to all modes)
        if options.get("palm_annot", True):
            final_files.append(
                #f"{config['output_dir']}/{sample}/duskmatter/{sample}_RdRp.tsv"
                f"{config['output_dir']}/{sample}/report/{sample}_Report_Diversity.html"
                #f"{config['output_dir']}{{sample}}/report/{{sample}}_Report_Diversity.html
            )
            
    return final_files

rule all:
    input: get_final_targets


# --- 6. Cleanup Handler (MODIFIED FOR MIXED MODES) ---

def cleanup_downloaded_fastqs(log):
    # Use the global flag to check if *any* SRA download happened
    download_only = config.get("options", {}).get("download_only", True)
    
    if ANY_SRA_DOWNLOADED and not download_only:
        print("\n--- INITIATING POST-SUCCESS CLEANUP ---")
        print("Workflow complete. Removing downloaded FASTQ files from 'raw_data/'...")
        
        files_removed = 0
        for sample, mode in SAMPLE_MODES.items():
            # ONLY clean up if the sample was processed in SRA mode
            if mode == 'SRA':
                f1 = os.path.join(RAW_DATA_DIR, f"{sample}_1.fastq.gz")
                f2 = os.path.join(RAW_DATA_DIR, f"{sample}_2.fastq.gz")
                
                if os.path.exists(f1):
                    os.remove(f1)
                    print(f"Removed: {f1}")
                    files_removed += 1
                if os.path.exists(f2):
                    os.remove(f2)
                    print(f"Removed: {f2}")
                    files_removed += 1
                
        print(f"Cleanup complete. {files_removed} files removed.")
        print("---------------------------------------\n")
        
    elif ANY_SRA_DOWNLOADED and download_only:
        print("\nWorkflow (download_only) complete. *Keeping* SRA downloaded FASTQ files.")
        
    else: # If no SRA was involved (only local FASTQ or CONTIGS)
        print("\nWorkflow complete. *Keeping* user-provided local FASTQ files.")

# Register the function to run on success
onsuccess:
    cleanup_downloaded_fastqs

# Optional: Add a handler for errors
onerror:
    print("\nWorkflow encountered an error. No files from 'raw_data/' were cleaned up.")
