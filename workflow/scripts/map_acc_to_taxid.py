# Em scripts/map_acc_to_taxid.py
import sys
import os
import csv
import gzip
import logging

# Configuração de logging: Snakemake redireciona o stderr para o arquivo .log
logging.basicConfig(
    filename=str(snakemake.log[0]),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# -----------------------------------------------------
# 1. Obter Entradas do Snakemake
# -----------------------------------------------------

try:
    hit_file = snakemake.input["hit_file"]
    output_file = snakemake.output["ids"]
    taxid_map_file = snakemake.params["taxid_map"]
    # Usamos o resultado da função Python de checagem do header
    header_present = snakemake.params["header_check"](hit_file)
except Exception as e:
    logging.error(f"Erro ao carregar parâmetros do Snakemake: {e}")
    sys.exit(1)

logging.info(f"Iniciando mapeamento para {hit_file}")
logging.info(f"Arquivo de TaxID Map: {taxid_map_file}")
logging.info(f"Header do arquivo de hits DIAMOND presente: {header_present}")


# -----------------------------------------------------
# 2. Construir o Mapa Accession -> TaxID
# -----------------------------------------------------
def build_taxid_map(taxid_path):
    """Lê o arquivo prot.accession2taxid.gz ou .gz e cria o dicionário de mapeamento."""
    taxid_map = {}
    
    # Usa gzip.open se o arquivo terminar com .gz, senão usa open()
    open_func = gzip.open if taxid_path.endswith('.gz') else open
    read_mode = 'rt' # 'rt' for reading text, crucial for gzip.open

    try:
        logging.info(f"Lendo mapa de TaxID com {open_func.__name__}...")
        with open_func(taxid_path, read_mode) as f:
            # Assumimos que o TaxID Map está no formato: gi\taccession.version\ttaxid\tprotein_name
            # Precisamos da coluna 2 (accession.version) e 3 (taxid).
            reader = csv.reader(f, delimiter='\t')
            
            # Pula o header se estiver presente (assumindo que o primeiro campo é 'accession')
            try:
                first_line = next(reader)
                if first_line and first_line[1] == 'accession.version':
                    logging.info("TaxID Map: Header detectado e pulado.")
                else:
                    # Se não for header, processamos essa linha
                    if len(first_line) >= 3:
                        taxid_map[first_line[1]] = first_line[2]
            except StopIteration:
                 # Arquivo vazio ou com apenas o header
                 logging.warning("TaxID Map: Arquivo vazio ou com apenas o header.")
                 return taxid_map
            except Exception as e:
                 # Erro ao ler a primeira linha (e.g. compressão corrompida)
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


# -----------------------------------------------------
# 3. Processar Hits DIAMOND e Anexar TaxID
# -----------------------------------------------------

try:
    with open(hit_file, 'r', newline='') as infile, \
         open(output_file, 'w', newline='') as outfile:
        
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        processed_hits = 0
        
        for i, row in enumerate(reader):
            # Se for a primeira linha (i=0)
            if i == 0:
                # Se o header estiver presente (header_present=1), adicionamos a coluna 'taxid'
                if header_present:
                    writer.writerow(row + ['taxid'])
                    continue # Próxima linha
                else:
                    # Se não houver header, escrevemos a linha atual com o TaxID.
                    pass # Continua para a lógica de dados abaixo

            # Lógica de processamento de dados (Subject ID é a coluna 2, índice 1)
            if len(row) >= 2:
                protein_id = row[1]
                
                # Procura no mapa. Se não encontrar, assume 'NOT_FOUND'.
                # Nota: A chave no TaxID Map é sempre accession.version (ex: AAF30685.1)
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