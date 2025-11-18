# Snakefile

import os
from glob import glob

# --- 1. Load Configuration ---
configfile: "config/workflow_config.yaml"

# --- 2. Discover Samples ---
RAW_DATA_DIR = "raw_data"
# Use glob to find all '_1.fastq.gz' files, which represent paired-end samples
local_fastqs = glob(os.path.join(RAW_DATA_DIR, "*_1.fastq.gz"))

# This flag tracks if data was downloaded (True) or provided locally (False)
DOWNLOADED_SRA = False

if local_fastqs:
    print("Found local FASTQ files in 'raw_data/'. Deriving sample names from them.")
    # Extract sample names from filenames, e.g., "raw_data/SRR12345_1.fastq.gz" -> "SRR12345"
    SAMPLES = [os.path.basename(f).replace("_1.fastq.gz", "") for f in local_fastqs]
    # Perform a quick sanity check to ensure read pairs exist
    for sample in SAMPLES:
        paired_file = os.path.join(RAW_DATA_DIR, f"{sample}_2.fastq.gz")
        if not os.path.exists(paired_file):
            raise FileNotFoundError(f"Error: Found {sample}_1.fastq.gz but its pair ({paired_file}) is missing.")
else:
    print(" No local FASTQ files found. Attempting to read SRA accessions for download.")
    try:
        sra_file_path = config["sra_file"]
        with open(sra_file_path, "r") as f:
            SAMPLES = [line.strip() for line in f if line.strip()]
        if not SAMPLES:
            raise ValueError("The SRA file is empty.")
        
        # --- (NEW) MODIFICATION 2: Update flag ---
        # If we get here, it means we are going to download the data
        DOWNLOADED_SRA = True

    except (KeyError, FileNotFoundError):
        # This error is raised if 'sra_file' is not in the config or the path is wrong
        raise ValueError(
            "CRITICAL: No local FASTQ files found in 'raw_data/' AND the 'sra_file' path "
            "is either missing from config.yaml or the file path is invalid."
        )

print(f"Workflow will process the following samples: {SAMPLES}")

# --- 3. Define Global Environment Paths ---
DOWNLOAD = workflow.source_path(config["envs"]["download"])
FASTP = workflow.source_path(config["envs"]["fastp"])
SPADES = workflow.source_path(config["envs"]["spades"])
DIAMOND = workflow.source_path(config["envs"]["diamond"])
TAXONKIT = workflow.source_path(config["envs"]["taxonkit"])
SCRIPTS = workflow.source_path(config["envs"]["scripts"])
PANDOC = workflow.source_path(config["envs"]["pandoc"])
PALM = workflow.source_path(config["envs"]["palm_annot"])


# --- 4. Load the Workflow Rules ---
include: "workflow/rules/get_data.smk"
include: "workflow/rules/fastp.smk"
include: "workflow/rules/spades.smk"
include: "workflow/rules/diamond.smk"
include: "workflow/rules/lineage.smk"
include: "workflow/rules/reporting.smk"
include: "workflow/rules/duskmatter.smk"

def get_final_targets(wildcards):
    final_files = []
    options = config.get("options", {})
    download_only = options.get("download_only", False)

    if download_only:
        print("Workflow mode: DOWNLOAD ONLY")
        final_files.extend(
            expand(
                [
                    "raw_data/{sample}_1.fastq.gz",
                    "raw_data/{sample}_2.fastq.gz"
                ],
                sample=SAMPLES
            )
        )
    else:
        print("Workflow mode: FULL ANALYSIS")
        
        # --- FIX 1 ---
        # The filename {sample} must be escaped as {{sample}}
        final_files.extend(
            expand(
                f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_with_lineage.tsv",
                sample=SAMPLES
            )
        )
        
        if options.get("diversity_reads", False):
            # This one was already correct in your file, leaving it as-is
            final_files.extend(
                expand(
                    f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_reads_hits_with_lineage.tsv",
                    sample=SAMPLES
                )
            )
            
        if options.get("palm_annot", True):
            final_files.extend(
                expand(
                    f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.tsv",
                    sample=SAMPLES
                )
            )
    return final_files

rule all:
    input: get_final_targets
# --- 6. Cleanup Handler ---

def cleanup_downloaded_fastqs(log):
    """
    This function is called by Snakemake ONLY IF the
    entire workflow finishes successfully (onsuccess).
    """
    # Get the 'download_only' value from the config
    download_only = config.get("options", {}).get("download_only", True)

    # We ONLY delete if:
    # 1. The data was downloaded from SRA (DOWNLOADED_SRA == True)
    # 2. We are NOT in "download_only" mode (because the full analysis was done)
    
    if DOWNLOADED_SRA and not download_only:
        print("\n--- INITIATING POST-SUCCESS CLEANUP ---")
        print("Workflow complete. Removing downloaded FASTQ files from 'raw_data/'...")
        
        files_removed = 0
        for sample in SAMPLES:
            # Recreate the paths to the files that were downloaded
            f1 = os.path.join(RAW_DATA_DIR, f"{sample}_1.fastq.gz")
            f2 = os.path.join(RAW_DATA_DIR, f"{sample}_2.fastq.gz")
            
            # Remove the files, if they exist
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
        
    elif DOWNLOADED_SRA and download_only:
        print("\nWorkflow (download_only) complete. *Keeping* downloaded FASTQ files.")
        
    else: # If local_fastqs was True
        print("\nWorkflow complete. *Keeping* user-provided local FASTQ files.")

# Register the function to run on success
onsuccess:
    cleanup_downloaded_fastqs

# Optional: Add a handler for errors
onerror:
    print("\nWorkflow encountered an error. No files from 'raw_data/' were cleaned up.")