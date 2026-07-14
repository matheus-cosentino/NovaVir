# Tutorial avançado de configuração

Esta página reúne as principais opções para ajustar o pré-processamento com fastp e a escolha da ferramenta de montagem de contigs no DiscoVir.

## 1. Onde editar as configurações

As opções principais estão no arquivo [config/config.yaml](../config/config.yaml).

Ao alterar esse arquivo, salve as mudanças e rode novamente o pipeline. O fluxo irá ler as novas definições automaticamente.

## 2. Ajustando o fastp

O fastp é usado para filtrar e limpar os reads antes da montagem e das análises taxonômicas.

### Parâmetros principais

No bloco `fastp` do arquivo [config/config.yaml](../config/config.yaml), você pode ajustar:

- `length_required`: comprimento mínimo do read após o trim.
- `qualified_quality_phred`: limiar mínimo de qualidade Phred.

Exemplo:

```yaml
fastp:
  length_required:
    - 50
  qualified_quality_phred:
    - 30
```

### Quando ajustar cada parâmetro

- Se seus dados são de alta qualidade e você quer preservar reads mais curtos, pode aumentar `length_required`.
- Se os dados são mais ruidosos, pode reduzir `qualified_quality_phred` para valores como 20 ou 25.
- Para dados long-read, em geral é mais seguro usar um limite de tamanho mais conservador, por exemplo `500`.

### Recomendação prática

- Para Illumina padrão: mantenha `length_required` em `50` e `qualified_quality_phred` em `30`.
- Para amostras com baixa qualidade: experimente `25` ou `20` para qualidade e `30` ou `40` para tamanho mínimo.

## 3. Escolhendo a ferramenta de montagem de contigs

A escolha do assembler é controlada no bloco `tool.denovo` em [config/config.yaml](../config/config.yaml).

Exemplo:

```yaml
tool:
  denovo:
    - 'spades'
    # - 'megahit'
    # - 'flye'
```

### Quando usar cada opção

- `spades`: melhor opção para dados Illumina metagenômicos e metatranscriptômicos. É o padrão atual e costuma performar melhor para recuperação de contigs virais.
- `megahit`: mais rápido e geralmente útil quando você prioriza velocidade em vez de máxima qualidade de montagem.
- `flye`: indicado para dados long-read (Nanopore/PacBio). Não é a melhor opção para reads curtos Illumina.
- `raven`: outra opção para long-read, geralmente mais simples e mais rápida, mas com trade-offs de qualidade.

## 4. Ajustando o SPAdes

Se você optar por SPAdes, há dois pontos importantes para ajustar:

### 4.1 Algoritmo

No bloco `spades.algorithm` você define o modo de montagem:

```yaml
spades:
  algorithm:
    - '--meta'
```

- `--meta`: recomendado para dados metagenômicos.
- `--isolate`: pode ser melhor para amostras isoladas ou dados de alta cobertura.

### 4.2 K-mer

O bloco `spades.kmer` define os tamanhos de k-mer usados na montagem:

```yaml
spades:
  kmer:
    - 'auto'
```

- `auto`: é a opção mais simples e geralmente boa para começar.
- valores explícitos como `21`, `33`, `55` ou combinações como `21_33_55` podem ser usados para testes mais específicos.

## 5. Ajustando a montagem para long-read

Se você estiver trabalhando com dados long-read, prefira `flye` ou `raven` em vez de SPAdes.

No arquivo [config/config.yaml](../config/config.yaml), você pode ajustar o tipo de leitura para o Flye:

```yaml
flye:
  type: "nano-corr"
```

Você pode trocar para outras opções, como:

- `pacbio-raw`
- `pacbio-corr`
- `pacbio-hifi`
- `nano-raw`
- `nano-corr`
- `nano-hq`

## 6. Ajustando parâmetros downstream

Além do fastp e do assembler, outros parâmetros podem influenciar bastante o resultado:

- `diamond.min_contig_len`: tamanho mínimo de contig para anotação com DIAMOND.
- `diamond.evalue`: limiar de e-value para alinhamentos.
- `kraken2.confidence`: sensibilidade para classificação taxonômica.
- `palm_annot.minscore`: limiar de score para detecção de RdRp candidates.

## 7. Fluxo recomendado para começar

1. Comece com os valores padrão do arquivo [config/config.yaml](../config/config.yaml).
2. Se a qualidade dos reads for ruim, ajuste fastp.
3. Se você tiver Illumina curto, teste `spades` primeiro.
4. Se você tiver long-read, use `flye` ou `raven`.
5. Compare os resultados e ajuste um parâmetro por vez.

## 8. Exemplo de execução após editar a configuração

```bash
bash DiscoVir.sh --input <DIR> --output <DIR> --profile local --diamond --darkmatter
```

Se você estiver rodando em cluster, troque o profile para `profile_slurm`.
