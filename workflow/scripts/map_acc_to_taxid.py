# scripts/map_acc_to_taxid.py

import sys
import os
import csv
import gzip
import logging
import subprocess 

# Objeto Snakemake
s = snakemake 

# --- 1. Configuração e Obtenção de Entradas ---

# Configuração de logging: Usa o log do Snakemake (s.log[0])
logging.basicConfig(
    filename=str(s.log[0]),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

try:
    hit_file = s.input["hit_file"]
    output_file = s.output["ids"]
    taxid_map_file = s.params["taxid_map"]
    
    # 1.1. Checagem de Header (Lógica Interna, independente de params)
    # Verifica se a primeira linha do arquivo de hits começa com 'qseqid'
    with open(hit_file, 'r') as f:
        first_line_hit = f.readline()
        # Se começar com qseqid, header presente (1). Senão, não (0).
        header_present = 1 if first_line_hit.startswith('qseqid') else 0
        
except Exception as e:
    logging.error(f"Erro ao carregar parâmetros ou checar header do hit file: {e}")
    sys.exit(1)

logging.info(f"Iniciando mapeamento para {hit_file} (Header presente: {header_present})")
logging.info(f"Arquivo de TaxID Map: {taxid_map_file}")

# --- 2. Função de Construção do Mapa Accession -> TaxID ---

def build_taxid_map(taxid_path):
    """Lê o arquivo TaxID Map (GZ ou plain) e cria o dicionário de mapeamento."""
    taxid_map = {}
    
    # Usa gzip.open se o arquivo terminar com .gz, senão usa open() (Robusto)
    open_func = gzip.open if taxid_path.endswith('.gz') else open
    read_mode = 'rt' # 'rt' para leitura de texto, crucial para gzip.open

    try:
        logging.info(f"Lendo mapa de TaxID com {open_func.__name__}...")
        with open_func(taxid_path, read_mode) as f:
            reader = csv.reader(f, delimiter='\t')
            
            # Lógica para pular o header ou processar a primeira linha
            try:
                first_line = next(reader)
                # Verifica se a primeira linha é o cabeçalho padrão (e.g., 'accession.version')
                if first_line and len(first_line) >= 3 and first_line[1] == 'accession.version':
                    logging.info("TaxID Map: Header detectado e pulado.")
                else:
                    # Se não for header, processamos essa linha
                    if len(first_line) >= 3:
                        taxid_map[first_line[1]] = first_line[2]
            except StopIteration:
                 logging.warning("TaxID Map: Arquivo vazio.")
                 return taxid_map
            except Exception as e:
                 logging.error(f"TaxID Map: Erro ao ler a primeira linha (possível corrupção): {e}")
                 sys.exit(1)

            # Processa o resto do arquivo
            for row in reader:
                if len(row) >= 3:
                    # Mapeia accession.version (col 2) -> taxid (col 3)
                    taxid_map[row[1]] = row[2]
        
        logging.info(f"Mapa de TaxID construído com {len(taxid_map)} entradas.")
        return taxid_map
        
    except FileNotFoundError:
        logging.error(f"Arquivo TaxID Map não encontrado em: {taxid_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Erro crítico ao ler TaxID Map: {e}")
        sys.exit(1)


taxid_map = build_taxid_map(taxid_map_file)

# --- 3. Processar Hits DIAMOND e Anexar TaxID ---

try:
    with open(hit_file, 'r', newline='') as infile, \
         open(output_file, 'w', newline='') as outfile:
        
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        processed_hits = 0
        
        for i, row in enumerate(reader):
            
            # Se a primeira linha (i=0) for o header do hit file, escrevemos o novo header e pulamos
            if i == 0 and header_present == 1:
                writer.writerow(row + ['taxid'])
                continue 

            # Lógica de processamento de dados (Subject ID é a coluna 2, índice 1)
            if len(row) >= 2:
                protein_id = row[1]
                
                # Procura no mapa, se não encontrar, usa 'NOT_FOUND'.
                taxid = taxid_map.get(protein_id, "NOT_FOUND")
                
                writer.writerow(row + [taxid])
                processed_hits += 1
            else:
                 logging.warning(f"Linha {i+1} do hit file ignorada: formato inválido.")

    logging.info(f"Processamento concluído. {processed_hits} hits mapeados.")
    logging.info("Script finalizado com sucesso.")

except Exception as e:
    logging.error(f"Erro crítico durante o processamento do hit file: {e}")
    sys.exit(1)

# Sucesso!