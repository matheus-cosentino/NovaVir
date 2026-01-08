# scripts/map_acc_to_taxid.py

import sys
import os
import csv
import gzip
import logging

# Objeto Snakemake
s = snakemake 

# --- 1. Configuração e Logging ---
logging.basicConfig(
    filename=str(s.log[0]),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

try:
    hit_file = s.input["hit_file"]
    output_file = s.output["ids"]
    
    taxid_map_file_raw = s.input["taxid_map"]
    # Garante que seja string mesmo se vier como lista
    taxid_map_file = str(taxid_map_file_raw[0]) if isinstance(taxid_map_file_raw, list) else str(taxid_map_file_raw)

except Exception as e:
    logging.error(f"Erro ao carregar parâmetros: {e}")
    sys.exit(1)

logging.info(f"Processando arquivo de hits: {hit_file}")
logging.info(f"Usando banco de dados de mapeamento: {taxid_map_file}")

# --- 2. Passo 1: Identificar quais proteínas precisamos procurar ---
# Isso economiza memória, pois não carregamos o banco inteiro, apenas o necessário.

needed_accessions = set()
header_present = 0

try:
    with open(hit_file, 'r') as f:
        # Detectar header
        first_line = f.readline()
        if first_line.startswith('qseqid'):
            header_present = 1
        
        # Voltar ao inicio
        f.seek(0)
        
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i == 0 and header_present:
                continue
            if len(row) >= 2:
                # Subject ID (Proteína) é a coluna 2 (índice 1)
                needed_accessions.add(row[1])

    logging.info(f"Identificadas {len(needed_accessions)} proteínas únicas para buscar no banco.")

except Exception as e:
    logging.error(f"Erro ao ler arquivo de hits para pré-filtragem: {e}")
    sys.exit(1)

# --- 3. Passo 2: Construir o Mapa Filtrado (Apenas o necessário) ---

taxid_map = {}

open_func = gzip.open if taxid_map_file.endswith('.gz') else open
read_mode = 'rt' # rt = read text

try:
    logging.info("Lendo arquivo gigante de mapeamento (modo streaming)...")
    
    with open_func(taxid_map_file, read_mode) as f:
        reader = csv.reader(f, delimiter='\t')
        
        count = 0
        found_count = 0
        
        for row in reader:
            count += 1
            if count % 10000000 == 0:
                logging.info(f"Linhas processadas do banco: {count/1000000} Milhões...")

            # O formato esperado é: accession.version (col 1), taxid (col 2) 
            # (No arquivo prot.accession2taxid.gz: accession (0), accession.version (1), taxid (2), gi (3))
            # Vamos checar a coluna 1 (accession.version)
            if len(row) >= 3:
                acc_ver = row[1]
                
                if acc_ver in needed_accessions:
                    taxid_map[acc_ver] = row[2]
                    found_count += 1
                    
                    # Opcional: Se já achamos todas, podemos parar cedo (break)
                    if found_count == len(needed_accessions):
                        logging.info("Todas as proteínas foram encontradas. Parando leitura do banco.")
                        break
                        
    logging.info(f"Mapeamento concluído. {found_count} taxids recuperados de {len(needed_accessions)} necessários.")

except Exception as e:
    logging.error(f"Erro ao ler arquivo de TaxID Map: {e}")
    sys.exit(1)

# --- 4. Passo 3: Escrever o arquivo final ---

try:
    with open(hit_file, 'r', newline='') as infile, \
         open(output_file, 'w', newline='') as outfile:
        
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        for i, row in enumerate(reader):
            
            # Header
            if i == 0 and header_present == 1:
                writer.writerow(row + ['taxid'])
                continue 

            if len(row) >= 2:
                protein_id = row[1]
                # Busca no dicionário filtrado
                taxid = taxid_map.get(protein_id, "NOT_FOUND")
                writer.writerow(row + [taxid])

    logging.info("Arquivo de saída gerado com sucesso.")

except Exception as e:
    logging.error(f"Erro ao escrever arquivo final: {e}")
    sys.exit(1)