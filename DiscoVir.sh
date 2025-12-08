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
blu="\033[1;34m"  
green="\033[1;32m" 
red="\033[1;31m"   
ylo="\033[1;33m"   
nc="\033[0m"       

# Nome do ambiente Conda (deve bater com workflow/envs/DiscoVir.yaml)
ENV_NAME="DiscoVir"
ENV_FILE="workflow/envs/DiscoVir.yaml"

#####################
# --- Functions --- #
#####################

# --- Spinner --- #
run_with_spinner() {
    local pid
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null
    local sequence="GGACTATCCTAAGTGTGACCGTGCTTTGCCGAGCATGATTAGGATGATTTCTGCCATGATACTTGGCTCTAAGCACACAACTTGCTGCACAAATAGTGATAGGTATTACAGATTGTGCAATGAGTTGGCACAAGTGCTCACTGAAGTTGTTTATTCCAATGGTGGTTTTTATTTTAAACCAGGAGGTACAACTTCAGGTGATGCAACTACAGCATATGCCAATTCTGTTTTCAACATATTCCAGGCTGTCAGTGCTAACATTAACCGTTTGCTCACTGTTGACAGTTATGCTATTCATAATGATTCTGTCAAGAGTTTGCAGAGGCAGTTGTATGACAATTGCTACCGTGCCACTTCTGTA"
    local seq_len=${#sequence}
    local width=15
    local pos=0
    
    # Hide cursor
    tput civis
    local cA="\033[1;32m"
    local cT="\033[1;31m"
    local cC="\033[1;34m"
    local cG="\033[1;33m"
    local nc="\033[0m"

    while kill -0 $pid 2>/dev/null; do
        local chunk=""
        for (( j=0; j<width; j++ )); do
            local idx=$(( (pos + j) % seq_len ))
            local base="${sequence:$idx:1}"
            case "$base" in
                A) chunk="${chunk}${cA}A${nc}" ;;
                T) chunk="${chunk}${cT}T${nc}" ;;
                C) chunk="${chunk}${cC}C${nc}" ;;
                G) chunk="${chunk}${cG}G${nc}" ;;
                *) chunk="${chunk}${base}" ;;
            esac
        done
        printf "\r\033[K [ ${chunk} ] Processing..."
        pos=$(( (pos + 1) % seq_len ))
        sleep 0.1
    done
    wait $pid
    local exit_code=$?
    tput cnorm
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
 ${green}
 DiscoVir${nc}: Viral Metagenomics & 'Dusk Matter' Discovery

 ${green}Author${nc}: MSc. Matheus Cosentino 
 ${green}Version${nc}: 12.2025

 ${ylo}Usage:${nc}
  bash DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 ${ylo}Required Arguments:${nc}
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved

 ${ylo}Optional Arguments:${nc}
   --sra <FILE>         Text file containing SRA Accession IDs for download.
   --jobs <INT>         Number of jobs (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)

 ${ylo}Database Overrides (Use external DBs):${nc}
   --diamond_db <FILE>  External Diamond database (.dmnd).
   --kraken2 <DIR>      External Kraken2 database directory (Must contain hash.k2d, opts.k2d, taxo.k2d)
   --taxdump <DIR>      Directory containing nodes.dmp and names.dmp
   --taxmap <FILE>      Path to prot.accession2taxid.gz

 ${ylo}Module Toggles (Enable/Disable Analysis):${nc}
   --assembly           Enable De Novo Assembly (Default: False)
   --kraken2            Enable Kraken2 Taxonomy (Default: False)
   --diamond            Enable Diamond Taxonomy (Default: False)
   --duskmatter         Enable Palm Annot / Dusk Matter (Default: False)
   --remove-download    Do NOT keep downloaded SRA files (Default: Keep)
   --skip-reads-kraken  Skip Reads Kraken2 ID (Default: Run)
   --skip-reads-diamond Skip Reads Diamond ID (Default: Run)
 
 ${ylo}Flags:${nc}
   -h, --help           Show this help message
   -v, --version        Show version
 "
}

version(){
    echo "DiscoVir v.12.2025"
}

###################################
# --- Environment Management --- #
###################################

manage_environment(){
    echo -e "${blu}[INFO]${nc} Checking Conda environment '${ylo}${ENV_NAME}${nc}'..."
    __conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else
        if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
            . "${HOME}/miniconda3/etc/profile.d/conda.sh"
        elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
            . "${HOME}/anaconda3/etc/profile.d/conda.sh"
        else
            echo -e "${red}[ERROR]${nc} Could not find conda.sh. Ensure Conda is installed."
            exit 1
        fi
    fi

    if [[ ! -f "$ENV_FILE" ]]; then
         if [[ -f "workflow/envs/DiscoVir.yaml" ]]; then
            ENV_FILE="workflow/envs/DiscoVir.yaml"
         else
            echo -e "${red}[ERROR]${nc} Environment file $ENV_FILE not found!"
            exit 1
         fi
    fi

    if conda env list | grep -q "^${ENV_NAME} "; then
        echo -e "       > Environment found: ${green}Yes${nc}"
        echo -e "${blu}[INFO]${nc} Updating environment from ${ylo}$ENV_FILE${nc}..."
        run_with_spinner conda env update --name "$ENV_NAME" --file "$ENV_FILE" --prune --quiet
    else 
        echo -e "       > Environment found: ${red}No${nc}"
        echo -e "${blu}[INFO]${nc} Creating environment from ${ylo}$ENV_FILE${nc}..."
        run_with_spinner conda env create --name "$ENV_NAME" --file "$ENV_FILE" --quiet
    fi

    echo -ne "${blu}[INFO]${nc} Activating environment... "
    conda activate "$ENV_NAME"
    if [[ $? -eq 0 ]]; then 
        echo -e "${green}Active${nc}"
    else 
        echo -e "${red}Failed${nc}"
        exit 1
    fi
}

###################################
# --- Resource Management --- #
###################################

setup_resources(){
    # Diretório base de recursos do workflow
    RES_DIR="$workdir/resources"

    # --- 1. Diamond DB ---
    # O config.yaml deve apontar para: "resources/diamond/database.dmnd"
    if [[ -n "$diamond_db" ]]; then
        if [[ -f "$diamond_db" ]]; then
            echo -e "${blu}[INFO]${nc} Linking external Diamond DB to resources/diamond/..."
            
            # Cria a pasta se não existir
            mkdir -p "$RES_DIR/diamond"
            
            # Remove link ou arquivo antigo para evitar conflito
            rm -f "$RES_DIR/diamond/database.dmnd"
            
            # Cria o link simbólico com um nome padronizado "database.dmnd"
            ln -sf "$diamond_db" "$RES_DIR/diamond/database.dmnd"
        else
            echo -e "${red}[ERROR]${nc} Diamond DB file not found: $diamond_db"
            exit 1
        fi
    fi

    # --- 2. Kraken2 DB ---
    # O config.yaml deve apontar para: "resources/kraken2/"
    if [[ -n "$kraken2" ]]; then
        # Remove a barra final se houver
        kraken2=${kraken2%/}

        if [[ -d "$kraken2" ]]; then
            echo -e "${blu}[INFO]${nc} Linking external Kraken2 DB to resources/kraken2/..."
            
            # Remove o diretório ou link 'kraken2' existente dentro de resources
            # Atenção: Isso substitui a pasta local pelo link para a externa
            rm -rf "$RES_DIR/kraken2"
            
            # Cria o link simbólico do diretório inteiro
            ln -sfn "$kraken2" "$RES_DIR/kraken2"
        else
            echo -e "${red}[ERROR]${nc} Kraken2 directory not found: $kraken2"
            exit 1
        fi
    fi

    # --- 3. Taxonomy (Nodes & Names & Map) ---
    # O config.yaml deve apontar para: "resources/taxonomy/..."
    if [[ -n "$taxdump" ]] || [[ -n "$taxmap" ]]; then
        mkdir -p "$RES_DIR/taxonomy"
    fi

    if [[ -n "$taxdump" ]]; then
        if [[ -d "$taxdump" ]]; then
            echo -e "${blu}[INFO]${nc} Linking Taxonomy Dump to resources/taxonomy/..."
            
            if [[ -f "$taxdump/nodes.dmp" ]]; then
                ln -sf "$taxdump/nodes.dmp" "$RES_DIR/taxonomy/nodes.dmp"
            else 
                echo -e "${ylo}[WARN]${nc} nodes.dmp not found in $taxdump"
            fi
            
            if [[ -f "$taxdump/names.dmp" ]]; then
                ln -sf "$taxdump/names.dmp" "$RES_DIR/taxonomy/names.dmp"
            else
                 echo -e "${ylo}[WARN]${nc} names.dmp not found in $taxdump"
            fi
        else
            echo -e "${red}[ERROR]${nc} Taxdump directory not found: $taxdump"
            exit 1
        fi
    fi

    if [[ -n "$taxmap" ]]; then
        if [[ -f "$taxmap" ]]; then
            echo -e "${blu}[INFO]${nc} Linking Accession Map to resources/taxonomy/..."
            ln -sf "$taxmap" "$RES_DIR/taxonomy/prot.accession2taxid.gz"
        else
             echo -e "${red}[ERROR]${nc} Taxmap file not found: $taxmap"
             exit 1
        fi
    fi
}


###################################
# --- Sample Generation --- #
###################################

generate_sample_list(){
    echo -ne "${blu}[INFO]${nc} Generating sample list..."
    sample_list="$output/samples_detected.txt"
    mkdir -p "$output"
    
    if [[ -d "$input" ]]; then
        find "$input" -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) \
        | sed 's|.*/||' \
        | sed -E 's/(_1|_2|_R1|_R2|_unp|_orphans)?\.(fastq\.gz|fastq|fasta|fa|fas)$//' \
        > "$sample_list"
    else
        touch "$sample_list"
    fi

    if [[ -n "$sra" && -f "$sra" ]]; then
        echo -e "\n${ylo}[INFO]${nc} Appending SRA Accessions..."
        echo "" >> "$sample_list" 
        cat "$sra" >> "$sample_list"
        echo "" >> "$sample_list"
    fi

    sort -u "$sample_list" | sed '/^$/d' > "${sample_list}.tmp" && mv "${sample_list}.tmp" "$sample_list"
    count=$(wc -l < "$sample_list")
    
    if [[ "$count" -eq 0 ]]; then
        echo -e " ${red}FAILED${nc}"
        echo -e "${red}Error: No samples found locally and no SRA file provided.${nc}"
        exit 1
    else
        echo -e " ${green}OK ($count unique samples)${nc}"
    fi
}

###########################
# --- Main Execution --- #
###########################

workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
set -e -o pipefail

# Defaults
input="data"
jobs=15
profile="profile_slurm"

# --- Module Defaults (must match config.yaml structure) ---
# We use Python-style boolean strings ("True"/"False")
mod_keep_download="True"
mod_assembly="False"
mod_kraken2="False"
mod_diamond="False"
mod_duskmatter="False"
mod_reads_kraken2="True"
mod_reads_diamond="True"

# --- Argument Parsing --- #
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --sra) sra="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        --profile) profile="$2"; shift 2 ;;
        
        ## Database Overrides
        --diamond_db) diamond_db="$2"; shift 2 ;;
        --kraken2) kraken2="$2"; shift 2 ;;
        --taxdump) taxdump="$2"; shift 2 ;;
        --taxmap) taxmap="$2"; shift 2 ;;

        # --- Module Toggles (Flags) ---
        --assembly) mod_assembly="True"; shift ;;
        --kraken2) mod_kraken2="True"; shift ;;
        --diamond) mod_diamond="True"; shift ;;
        --duskmatter) mod_duskmatter="True"; shift ;;
        
        # --- Negative Flags (Disable defaults) ---
        --remove-download) mod_keep_download="False"; shift ;;
        --skip-reads-kraken) mod_reads_kraken2="False"; shift ;;
        --skip-reads-diamond) mod_reads_diamond="False"; shift ;;

        # --- Others Flafs --- #
        -h|--help) help; exit 0 ;;
        -v|--version) version; exit 0 ;;
        *) echo -e "${red}[ERROR]${nc} Unknown argument: $1"; help; exit 1 ;;
    esac
done

# Basic Validation
if [[ -z "$output" ]]; then 
    echo -e "${red}[ERROR]${nc} Output directory is required. Use --output <DIR>"
    exit 1
fi

manage_environment
setup_resources
generate_sample_list

echo -e "${blu}[INFO]${nc} Initializing DiscoVir Workflow..."

# Config APENAS para inputs de dados, sem sobrescrever resources
modules_config="{'keep_download':${mod_keep_download},'assembly':${mod_assembly},'kraken2':${mod_kraken2},'diamond':${mod_diamond},'duskmatter':${mod_duskmatter},'reads_kraken2':${mod_reads_kraken2},'reads_diamond':${mod_reads_diamond}}"
config_override="data_dir='$input' output_dir='$output' sample_list='$sample_list' modules=$modules_config"

# --- Workflow Execution Steps ---

echo -e "\n${green}> Snakemake: Unlocking working directory...${nc}"
snakemake --profile $profile --config $config_override --unlock --quiet

echo -e "\n${green}> Snakemake: Creating conda environments (if needed)...${nc}"
snakemake --profile $profile --config $config_override --use-conda --conda-create-envs-only --quiet

echo -e "\n${green}> Snakemake: Performing a dry-run...${nc}"
snakemake --profile $profile --jobs $jobs --use-conda --config $config_override --dry-run
echo "---------------------------------------------------"

echo -e "\n${green}> Snakemake: Starting main execution...${nc}"
snakemake --profile $profile \
    --jobs $jobs \
    --use-conda \
    --config $config_override \
    --keep-going 

echo -e "\n${blu}[INFO]${nc} Deactivating ${ylo}${ENV_NAME}${nc} conda environment."
conda deactivate
echo -e "${green}✔ Done.${nc}"