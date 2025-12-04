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
    local spinner=("A" "T" "G" "T" "G" "T" "T" "C" "T" "G" "A" "C" "A" "A" "C" "A" "C" "G" "A" "T" "C" "A" "A" "C" "A" "T" "G")    local i=0
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
 ${green}
 DiscoVir${nc}: Viral Metagenomics & 'Dusk Matter' Discovery

 ${green}Author${nc}: MSc. Matheus Cosentino 
 
 ${green}Version${nc}: 12.2025

 ${ylo}
 Usage: ${nc}
  bash DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 ${ylo} Required Arguments:${nc}
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved

 ${ylo} Optional Arguments:${nc}
   --sra <FILE>         Text file containing SRA Accession IDs for download.
   --jobs <INT>         Number of jobs (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)

 ${ylo} Database Overrides (Use external DBs):${nc}
   --diamond_db <FILE>  External Diamond database (.dmnd).
                        (Script auto-detects taxonomy files in the same directory)
   --kraken2 <DIR>      External Kraken2 database directory
 
 ${ylo} Flags:${nc}
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

    # --- FIX: Always initialize conda for this script session ---
    # We try to locate the conda setup in standard locations and source it.
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

    # Check if env exists
    if conda env list | grep -q "^${ENV_NAME} "; then
        echo -e "       > Environment found: ${green}Yes${nc}"
    else 
        echo -e "       > Environment found: ${red}No${nc} (Will attempt to create it via Snakemake)"
    fi

    # Activate
    echo -ne "${blu}[INFO]${nc} Activating environment... "
    conda activate "$ENV_NAME"
    
    if [[ $? -eq 0 ]]; then 
        echo -e "${green}Active${nc}"
    else 
        echo -e "${red}Failed${nc}"
        echo -e "${ylo}Tip: If the environment doesn't exist yet, Snakemake will handle it.${nc}"
        # We don't exit here because Snakemake might create it later
    fi
}

###################################
# --- Resource Management --- #
###################################

setup_resources(){
    LINK_DIR="$workdir/resources/links"
    mkdir -p "$LINK_DIR"
    
    # 1. Diamond DB & Taxonomy Auto-Detection
    if [[ -n "$diamond_db" ]]; then
        if [[ -f "$diamond_db" ]]; then
            echo -e "${blue}[INFO]${nc} Linking external Diamond DB..."
            ln -sf "$diamond_db" "$LINK_DIR/external_diamond.dmnd"
            
            # Variável para a chave 'diamond'
            res_diamond="'$LINK_DIR/external_diamond.dmnd'"

            # AUTO-DETECTION: Procura arquivos de taxonomia na mesma pasta do .dmnd
            DB_DIR=$(dirname "$diamond_db")
            
            # Para chave 'taxonnodes'
            if [[ -f "$DB_DIR/nodes.dmp" ]]; then
                echo -e "       > Auto-detected ${green}nodes.dmp${nc}"
                ln -sf "$DB_DIR/nodes.dmp" "$LINK_DIR/nodes.dmp"
                res_nodes="'$LINK_DIR/nodes.dmp'"
            fi
            
            # Para chave 'taxonnames'
            if [[ -f "$DB_DIR/names.dmp" ]]; then
                echo -e "       > Auto-detected ${green}names.dmp${nc}"
                ln -sf "$DB_DIR/names.dmp" "$LINK_DIR/names.dmp"
                res_names="'$LINK_DIR/names.dmp'"
            fi
            
            # Para chave 'taxonmap'
            if [[ -f "$DB_DIR/prot.accession2taxid.gz" ]]; then
                echo -e "       > Auto-detected ${green}prot.accession2taxid.gz${nc}"
                ln -sf "$DB_DIR/prot.accession2taxid.gz" "$LINK_DIR/prot.accession2taxid.gz"
                res_map="'$LINK_DIR/prot.accession2taxid.gz'"
            fi

        else
            echo -e "${red}[WARNING]${nc} Diamond DB not found: $diamond_db. Using default."
        fi
    fi

    # 2. Kraken2 DB Override
    if [[ -n "$kraken2" ]]; then
        if [[ -d "$kraken2" ]]; then
             echo -e "${blue}[INFO]${nc} Linking external Kraken2 DB..."
            rm -f "$LINK_DIR/external_kraken2"
            ln -sfn "$kraken2" "$LINK_DIR/external_kraken2"
            # Variável para a chave 'kraken2'
            res_kraken="'$LINK_DIR/external_kraken2'"
        else
            echo -e "${red}[WARNING]${nc} Kraken2 dir not found: $kraken2. Using default."
        fi
    fi
}

###################################
# --- Sample Generation --- #
###################################

generate_sample_list(){
    echo -ne "${blue}[INFO]${nc} Generating sample list..."
    sample_list="$output/samples_detected.txt"
    mkdir -p "$output"
    
    # 1. Local Files
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
        echo -e "\n${ylo}[INFO]${nc} Appending SRA Accessions..."
        cat "$sra" >> "$sample_list"
        echo "" >> "$sample_list"
    fi

    sed -i '/^$/d' "$sample_list"
    count=$(wc -l < "$sample_list")
    
    if [[ "$count" -eq 0 ]]; then
        echo -e " ${red}FAILED${nc}"
        echo -e "${red}Error: No samples found locally and no SRA file provided.${nc}"
        exit 1
    else
        echo -e " ${green}OK ($count samples)${nc}"
    fi
}

###########################
# --- Main Execution --- #
###########################

workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
set -e -o pipefail

# Parse Arguments
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
   shift
done

if [[ -n "$help" || -n "$h" ]]; then help; exit 0; fi
if [[ -n "$version" || -n "$v" ]]; then version; exit 0; fi

if [[ -z "$output" ]]; then help; exit 1; fi
input=${input:-"data"} 
profile=${profile:-profile_slurm}
jobs=${jobs:-15}

manage_environment
setup_resources
generate_sample_list

echo -e "${blue}[INFO]${nc} Initializing DiscoVir Workflow..."

# Configuração Base (Caminhos e Lista de Amostras)
config_override="data_dir='$input' output_dir='$output' sample_list='$sample_list'"

# Construção do Override de Recursos (Mapeamento exato para o config.yaml)
# O Snakemake aceita: resources={key1:val1, key2:val2}
resource_str=""

# Mapeia variáveis do bash para as chaves do config.yaml
if [[ -n "$res_diamond" ]]; then resource_str="${resource_str}diamond:$res_diamond,"; fi
if [[ -n "$res_kraken" ]];  then resource_str="${resource_str}kraken2:$res_kraken,"; fi
if [[ -n "$res_nodes" ]];   then resource_str="${resource_str}taxonnodes:$res_nodes,"; fi
if [[ -n "$res_names" ]];   then resource_str="${resource_str}taxonnames:$res_names,"; fi
if [[ -n "$res_map" ]];     then resource_str="${resource_str}taxonmap:$res_map,"; fi

# Remove a última vírgula
resource_str=${resource_str%,}

# Se montamos alguma string de recurso, adicionamos ao comando
if [[ -n "$resource_str" ]]; then
    config_override="$config_override resources={$resource_str}"
fi

# --- Workflow Execution Steps ---

# 1. Unlock Directory (in case of previous errors)
echo -e "\n${green}> Snakemake: Unlocking working directory...${nc}"
snakemake --profile $profile --config $config_override --unlock --quiet

# 2. Create Conda Environments
echo -e "\n${green}> Snakemake: Creating conda environments (if needed)...${nc}"
snakemake --profile $profile --config $config_override --use-conda --conda-create-envs-only --quiet

# 3. Dry-Run
echo -e "\n${green}> Snakemake: Performing a dry-run...${nc}"
snakemake --profile $profile --jobs $jobs --use-conda --config $config_override --dry-run
echo "---------------------------------------------------"

# 4. Final Execution
echo -e "\n${green}> Snakemake: Starting main execution...${nc}"
snakemake --profile $profile \
    --jobs $jobs \
    --use-conda \
    --config $config_override \
    --keep-going \
    --rerun-incomplete

###############################################################################
### DEACTIVATE WORKFLOW-CORE ###
#################################

# Deactivate workflow-core conda environment.
echo -e "\n${blue}[INFO]${nc} Deactivating ${ylo}${ENV_NAME}${nc} conda environment."
conda deactivate
echo -e "${green}✔ Done.${nc}"