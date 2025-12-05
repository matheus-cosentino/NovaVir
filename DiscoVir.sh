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
    # Run the command in the background
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null

    # Define colors for bases (A=Green, T=Red, C=Blue, G=Yellow)
    local cA="\033[1;32m" # Green
    local cT="\033[1;31m" # Red
    local cC="\033[1;34m" # Blue
    local cG="\033[1;33m" # Yellow
    local nc="\033[0m"    # No Color

    # The full DNA sequence to scroll
    local sequence="GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAATCTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGTTCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAAACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTGCCAGAACACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATATCCCTCTAGAGGAGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATACCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAGACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTGAGAGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTCAAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGTCGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTCAACTTGACCGCTCTGAGCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAGTTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTCCAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTAACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCACTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATCACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACTATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTCTT"
    local seq_len=${#sequence}
    
    # Window size (how many base pairs to show at once)
    local width=15
    local pos=0

    # Hide cursor
    tput civis

    while kill -0 $pid 2>/dev/null; do
        # Construct the visible window
        local chunk=""
        for (( j=0; j<width; j++ )); do
            # Calculate index with wrapping
            local idx=$(( (pos + j) % seq_len ))
            local base="${sequence:$idx:1}"
            
            # Colorize the base
            case "$base" in
                A) chunk="${chunk}${cA}A${nc}" ;;
                T) chunk="${chunk}${cT}T${nc}" ;;
                C) chunk="${chunk}${cC}C${nc}" ;;
                G) chunk="${chunk}${cG}G${nc}" ;;
                *) chunk="${chunk}${base}" ;;
            esac
        done

        # Print the sliding window
        printf "\r\033[K [ ${chunk} ] Processing..."
        
        # Advance position
        pos=$(( (pos + 1) % seq_len ))
        
        # Speed of animation
        sleep 0.1
    done
    
    wait $pid
    local exit_code=$?
    
    # Restore cursor
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
   --taxdump <DIR>      Directory containing nodes.dmp and names.dmp
   --taxmap <FILE>      Path to prot.accession2taxid.gz
 
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

    # 1. Initialize Conda for this script session
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

    # 2. Verify YAML file exists
    if [[ ! -f "$ENV_FILE" ]]; then
         if [[ -f "workflow/envs/DiscoVir.yaml" ]]; then
            ENV_FILE="workflow/envs/DiscoVir.yaml"
         else
            echo -e "${red}[ERROR]${nc} Environment file $ENV_FILE not found!"
            exit 1
         fi
    fi

    # 3. Logic: Check -> Update (with Spinner) -> Activate
    if conda env list | grep -q "^${ENV_NAME} "; then
        echo -e "       > Environment found: ${green}Yes${nc}"
        
        echo -e "${blu}[INFO]${nc} Updating environment from ${ylo}$ENV_FILE${nc}..."
        
        # --- SPINNER APPLIED HERE ---
        # The spinner function runs the command in the background and handles exit codes
        run_with_spinner conda env update --name "$ENV_NAME" --file "$ENV_FILE" --prune --quiet
        
    else 
        echo -e "       > Environment found: ${red}No${nc}"
        echo -e "${blu}[INFO]${nc} Creating environment from ${ylo}$ENV_FILE${nc}..."
        
        # --- SPINNER APPLIED HERE TOO (Recommended) ---
        run_with_spinner conda env create --name "$ENV_NAME" --file "$ENV_FILE" --quiet
    fi

    # 4. Final Activation
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
    LINK_DIR="$workdir/resources/links"
    mkdir -p "$LINK_DIR"
    
    # --- 1. Diamond DB & Auto-Detection ---
    if [[ -n "$diamond_db" ]]; then
        if [[ -f "$diamond_db" ]]; then
            echo -e "${blu}[INFO]${nc} Linking external Diamond DB..."
            ln -sf "$diamond_db" "$LINK_DIR/external_diamond.dmnd"
            res_diamond="'$LINK_DIR/external_diamond.dmnd'"

            # AUTO-DETECTION (Lower priority)
            # Only runs if explicit tax flags are NOT provided later
            DB_DIR=$(dirname "$diamond_db")
            
            if [[ -f "$DB_DIR/nodes.dmp" ]]; then
                ln -sf "$DB_DIR/nodes.dmp" "$LINK_DIR/nodes.dmp"
                res_nodes="'$LINK_DIR/nodes.dmp'"
            fi
            
            if [[ -f "$DB_DIR/names.dmp" ]]; then
                ln -sf "$DB_DIR/names.dmp" "$LINK_DIR/names.dmp"
                res_names="'$LINK_DIR/names.dmp'"
            fi
            
            if [[ -f "$DB_DIR/prot.accession2taxid.gz" ]]; then
                ln -sf "$DB_DIR/prot.accession2taxid.gz" "$LINK_DIR/prot.accession2taxid.gz"
                res_map="'$LINK_DIR/prot.accession2taxid.gz'"
            fi
        else
            echo -e "${red}[WARNING]${nc} Diamond DB not found: $diamond_db. Using default."
        fi
    fi

    # --- 2. Kraken2 DB Override ---
    if [[ -n "$kraken2" ]]; then
        if [[ -d "$kraken2" ]]; then
             echo -e "${blu}[INFO]${nc} Linking external Kraken2 DB..."
            rm -f "$LINK_DIR/external_kraken2"
            ln -sfn "$kraken2" "$LINK_DIR/external_kraken2"
            res_kraken="'$LINK_DIR/external_kraken2'"
        else
            echo -e "${red}[WARNING]${nc} Kraken2 dir not found: $kraken2. Using default."
        fi
    fi

    # --- 3. TaxDump Override (High Priority) ---
    # This handles nodes.dmp and names.dmp
    if [[ -n "$taxdump" ]]; then
        if [[ -d "$taxdump" ]]; then
            echo -e "${blu}[INFO]${nc} Linking Taxonomy Dump from: ${ylo}$taxdump${nc}"
            
            # Link nodes.dmp
            if [[ -f "$taxdump/nodes.dmp" ]]; then
                ln -sf "$taxdump/nodes.dmp" "$LINK_DIR/nodes.dmp"
                res_nodes="'$LINK_DIR/nodes.dmp'"
                echo -e "       > Linked ${green}nodes.dmp${nc}"
            else
                echo -e "${red}[ERROR]${nc} nodes.dmp not found in $taxdump"
                exit 1
            fi

            # Link names.dmp
            if [[ -f "$taxdump/names.dmp" ]]; then
                ln -sf "$taxdump/names.dmp" "$LINK_DIR/names.dmp"
                res_names="'$LINK_DIR/names.dmp'"
                echo -e "       > Linked ${green}names.dmp${nc}"
            else
                echo -e "${red}[ERROR]${nc} names.dmp not found in $taxdump"
                exit 1
            fi
        else
            echo -e "${red}[ERROR]${nc} Taxdump directory not found: $taxdump"
            exit 1
        fi
    fi

    # --- 4. TaxMap Override (High Priority) ---
    # This handles prot.accession2taxid.gz
    if [[ -n "$taxmap" ]]; then
        if [[ -f "$taxmap" ]]; then
            echo -e "${blu}[INFO]${nc} Linking Accession Map from: ${ylo}$taxmap${nc}"
            ln -sf "$taxmap" "$LINK_DIR/prot.accession2taxid.gz"
            res_map="'$LINK_DIR/prot.accession2taxid.gz'"
        else
            echo -e "${red}[ERROR]${nc} Accession map file not found: $taxmap"
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
    
    # 1. Arquivos Locais (Sobrescreve/Reseta o arquivo)
    if [[ -d "$input" ]]; then
        find "$input" -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) \
        | sed 's|.*/||' \
        | sed -E 's/(_1|_2|_R1|_R2|_unp|_orphans)?\.(fastq\.gz|fastq|fasta|fa|fas)$//' \
        > "$sample_list"
    else
        touch "$sample_list"
    fi

    # 2. Injeção de SRA (Com tratamento de quebra de linha)
    if [[ -n "$sra" && -f "$sra" ]]; then
        echo -e "\n${ylo}[INFO]${nc} Appending SRA Accessions..."
        
        # Garante quebra de linha antes de adicionar
        echo "" >> "$sample_list" 
        
        # Adiciona o conteúdo do SRA
        cat "$sra" >> "$sample_list"
        
        # Garante quebra de linha DEPOIS de adicionar (caso o txt do usuário não tenha)
        echo "" >> "$sample_list"
    fi

    # 3. Limpeza Final (Remove linhas vazias e duplicatas)
    # O comando 'sort -u' remove duplicatas e organiza a lista
    sort -u "$sample_list" | sed '/^$/d' > "${sample_list}.tmp" && mv "${sample_list}.tmp" "$sample_list"

    # Validação
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
    --keep-going 

###############################################################################
### DEACTIVATE WORKFLOW-CORE ###
#################################

# Deactivate workflow-core conda environment.
echo -e "\n${blue}[INFO]${nc} Deactivating ${ylo}${ENV_NAME}${nc} conda environment."
conda deactivate
echo -e "${green}✔ Done.${nc}"