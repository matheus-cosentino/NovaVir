import pandas as pd
import sqlite3
import re
import urllib.request
import json
import time
import sys

# Snakemake Configurações
ARQUIVO_TBL = snakemake.input.tbl
BANCO_SQLITE = snakemake.input.sqlite
SAIDA_CSV = snakemake.output.csv
LOG_FILE = snakemake.log[0]

sys.stdout = open(LOG_FILE, "w")
sys.stderr = sys.stdout

# Palavras que queremos ignorar para obter uma descrição melhor
BLACKLIST = {'protein', 'hypothetical', 'unnamed', 'predicted', 'unknown', 'virus', 'sequence', 'structure', 'homolog', 'orf'}

# Cache para não perguntar ao NCBI a mesma coisa várias vezes (mais rápido)
tax_cache = {}

def get_scientific_name(taxid):
    """Consulta a API do NCBI para traduzir o TaxID em Nome Científico"""
    if taxid == "N/A" or not taxid: return "N/A"
    if taxid in tax_cache: return tax_cache[taxid]
    
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={taxid}&retmode=json"
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            name = data['result'][str(taxid)]['scientificname']
            tax_cache[taxid] = name
            time.sleep(0.1) # Evita ser bloqueado pelo NCBI
            return name
    except:
        return f"TaxID:{taxid}"

def extrair_id_familia(nome_alvo):
    match = re.search(r'FAM(\d+)', str(nome_alvo))
    return int(match.group(1)) if match else None

def buscar_metadados(fam_id, conn):
    # 1. Buscar keywords e filtrar as genéricas
    query_kw = f"""
        SELECT k.str FROM keyword k
        JOIN fam_kw_seqnames fks ON k.id = fks.kwId
        WHERE fks.famId = {fam_id}
        ORDER BY fks.freq DESC
    """
    kws = pd.read_sql_query(query_kw, conn)
    
    descricao = "N/A"
    for word in kws['str']:
        if word.lower() not in BLACKLIST:
            descricao = word
            break
            
    # 2. Buscar o TaxID (usando o LCA da tabela family para maior precisão taxonómica)
    query_tax = f"SELECT LCAtaxid FROM family WHERE id = {fam_id}"
    tax = pd.read_sql_query(query_tax, conn)
    taxid = tax['LCAtaxid'].iloc[0] if not tax.empty else None
    
    return descricao, taxid

print(f"-> A processar resultados de: {ARQUIVO_TBL}...")
cols_hmmer = ["target_name", "query_name", "e_value", "score"]
df = pd.read_csv(ARQUIVO_TBL, comment='#', header=None, sep=r'\s+', usecols=[0,2,4,5], names=cols_hmmer)

conn = sqlite3.connect(BANCO_SQLITE)
results = []

for idx, row in df.iterrows():
    fam_id = extrair_id_familia(row['target_name'])
    desc, tid = buscar_metadados(fam_id, conn)
    nome_sci = get_scientific_name(tid)
    
    results.append({
        'Query': row['query_name'],
        'Family': row['target_name'],
        'E-value': row['e_value'],
        'Description': desc,
        'TaxID': tid,
        'Scientific_Name': nome_sci
    })
    if idx % 10 == 0 and idx > 0: 
        print(f"Processados {idx} hits...")

df_final = pd.DataFrame(results)
df_final.to_csv(SAIDA_CSV, index=False)
conn.close()

print(f"-> Concluído! Verifique o ficheiro: {SAIDA_CSV}")
sys.stdout.close()
