#!/bin/bash
#SBATCH --job-name=DiscoVir
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs_slurm/DiscoVir_%j.out

echo "
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
##                                                      ##
##      ▗▄▄▄ ▗▄▄▄▖ ▗▄▄▖ ▗▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄        ##
##      ▐▌  █  █  ▐▌   ▐▌   ▐▌ ▐▌▐▌  ▐▌  █  ▐▌ ▐▌       ##
##      ▐▌  █  █   ▝▀▚▖▐▌   ▐▌ ▐▌▐▌  ▐▌  █  ▐▛▀▚▖       ##
##      ▐▙▄▄▀▗▄█▄▖▗▄▄▞▘▝▚▄▄▖▝▚▄▞▘ ▝▚▞▘ ▗▄█▄▖▐▌ ▐▌       ##
##                                                      ##
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
## Name ___________________ DiscoVir.slurm              ##
## Version ________________ v.2025.12                   ##
## Author _________________ Matheus Cosentino           ##
## Affiliation ____________ IRD TransVIHMI & UFRJ PGGEN ##
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
"

# --- 1. Create a unique temporary directory on the node's local scratch ---
# $SLURM_JOB_ID is a unique ID given to your job by Slurm
# (Activate your main Conda environment that has Snakemake)
# Source the conda.sh script from your base environment to initialize conda
source /usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh


# --- 2. Run Snakemake using the shadow directory ---
# This tells Snakemake to run everything inside our unique temp directory
#bash DiscoVir.sh --input data --output results --diamond_db /share/banks/diamond/nr/21-09/diamond_nr.dmnd --taxdump /share/banks/taxdump --taxmap /databases/banks/ncbi/taxonomy/accession2taxid/prot.accession2taxid.FULL  --skip-reads-kraken --sra accessions.txt
 
bash DiscoVir.sh --input data --output results  --skip-reads-kraken --skip-reads-diamond --assembly --diamond --taxmap /databases/banks/ncbi/taxonomy/accession2taxid/prot.accession2taxid.FULL  --diamond_db /share/banks/diamond/nr/21-09/diamond_nr.dmnd
