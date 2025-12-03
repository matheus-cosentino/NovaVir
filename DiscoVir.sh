#!/bin/bash
###################################################################################
#                               DiscoVir.sh                                       #
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

# --- Color Palettes ---
blue="\033[1;34m"  
green="\033[1;32m" 
red="\033[1;31m"   
ylo="\033[1;33m"   
nc="\033[0m"       

#####################
# --- Functions --- #
#####################

# --- Spinner --- #
run_with_spinner() {
    local pid
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null
    local spinner=("A" "T" "G" "C")
    local i=0
    while kill -0 $pid 2>/dev/null; do
        printf "\r\033[K [${green}${spinner[i]}${nc}] Processing...    "
        i=$(( (i+1) % 4 ))
        sleep 0.2
    done
    wait $pid
    local exit_code=$?
    printf "\r\033[K"
    if [ $exit_code -eq 0 ]; then
        echo -e "${green}✔ Job done!${nc}"
    else
        echo -e "${red}✖ Job failed with exit code $exit_code.${nc}"
        exit $exit_code
    fi
}

# --- Help --- #
help(){
 echo -e "
 ${blue}DiscoVir: Viral Metagenomics & 'Dusk Matter' Discovery${nc}
 Author: MSc. Matheus Cosentino | Version: 12.2025 

 Usage: 
   srun DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 Required Arguments:
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved

 Optional Arguments:
   --sra <FILE>         Text file containing SRA Accession IDs (one per line) for download.
   --jobs <INT>         Number of jobs (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)

 Database Overrides:
   --diamond_db <FILE>  Path to Diamond database (.dmnd)
   --kraken2 <DIR>      Directory containing Kraken2 database
 
 Flags:
   -h, --help           Show this help message
   -v, --version        Show version
 "
}

version(){
    echo "DiscoVir v.2025.12"
}

###############################
# --- Check Prerequisites --- #
###############################

check_conda(){
    echo -ne "${blue}[INFO]${nc} Checking Conda availability..."
    if command -v conda &> /dev/null; then
        echo -e " ${green}OK${nc}"
    else
        echo -e " ${red}FAILED${nc}"
        echo "Error: Conda is not installed or not in PATH."
        exit 1
    fi
}

###################################
# --- Resource Management --- #
###################################

setup_resources(){
    echo -e "${blue}[INFO]${nc} Setting up resources symlinks..."
    
    mkdir -p "$workdir/resources_links"

    # 1. Diamond DB
    if [[ -n "$diamond_db" && -f "$diamond_db" ]]; then
        ln -sf "$diamond_db" "$workdir/resources_links/diamond.dmnd"
        diamond_path="$workdir/resources_links/diamond.dmnd"
    fi

    # 2. Kraken2 DB
    if [[ -n "$kraken2" && -d "$kraken2" ]]; then
        ln -sfn "$kraken2" "$workdir/resources_links/kraken2_db"
        kraken2_path="$workdir/resources_links/kraken2_db"
    fi
}

###################################
# --- Sample Generation --- #
###################################

generate_sample_list(){
    echo -ne "${blue}[INFO]${nc} Generating sample list..."
    
    # Final sample file
    sample_list="$output/samples_detected.txt"
    mkdir -p "$output"
    
    # 1. Local Files Detection
    if [[ -d "$input" ]]; then
        find "$input" -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) \
        | sed 's|.*/||' \
        | sed -E 's/(_1|_2|_R1|_R2|_unp|_orphans)?\.(fastq\.gz|fastq|fasta|fa|fas)$//' \
        | sort | uniq > "$sample_list"
    else
        touch "$sample_list"
    fi

    # 2. SRA Injection
    if [[ -n "$sra" && -f "$sra" ]]; then
        echo -e "\n${ylo}[INFO]${nc} Appending SRA Accessions from $sra..."
        cat "$sra" >> "$sample_list"
    elif [[ -n "$sra" ]]; then
        echo -e "\n${red}[WARNING]${nc} SRA file specified but not found: $sra"
    fi

    # 3. Final Validation
    sed -i '/^$/d' "$sample_list"
    
    count=$(wc -l < "$sample_list")
    
    if [[ "$count" -eq 0 ]]; then
        echo -e " ${red}FAILED${nc}"
        echo -e "${red}Error: No samples found locally in '$input' and no valid SRA file provided.${nc}"
        exit 1
    else
        echo -e " ${green}OK ($count samples queued)${nc}"
    fi
}

###########################
# --- Main Execution --- #
###########################

workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)

# Parse Arguments
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
   shift
done

if [[ -z "$output" ]]; then
    help
    exit 1
fi

# Input assumes 'data' if not provided
input=${input:-"data"} 

profile=${profile:-profile_slurm}
jobs=${jobs:-15}

check_conda
setup_resources
generate_sample_list

echo -e "${blue}[INFO]${nc} Initializing DiscoVir Workflow..."

# Base config override
config_args="data_dir=$input output_dir=$output sample_list=$sample_list"

# Resource overrides
if [[ -n "$diamond_path" ]]; then 
    config_args="$config_args resources={diamond:\"$diamond_path\"}" 
fi
if [[ -n "$kraken2_path" ]]; then 
    config_args="$config_args resources={kraken2:\"$kraken2_path\"}" 
fi

# Build command
cmd="snakemake --profile $profile \
    --jobs $jobs \
    --use-conda \
    --config $config_args"

echo -e "${ylo}Running command:${nc}"
echo "$cmd"
echo "---------------------------------------------------"

eval $cmd