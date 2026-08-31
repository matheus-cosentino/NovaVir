# Advanced configuration tutorial

This page gathers the main options for tuning preprocessing with fastp and choosing the contig assembler in NovaVir.

## 1. Where to edit configuration

The main settings are in the file [config/config.yaml](../config/config.yaml).

After you edit this file, save it and run the pipeline again. The workflow will automatically read the updated values.

## 2. Tuning fastp

fastp is used to filter and clean reads before assembly and taxonomic analysis.

### Main parameters

In the `fastp` section of [config/config.yaml](../config/config.yaml), you can adjust:

- `length_required`: minimum read length after trimming.
- `qualified_quality_phred`: minimum Phred quality threshold.

Example:

```yaml
fastp:
  length_required:
    - 50
  qualified_quality_phred:
    - 30
```

### When to change each parameter

- If your data are high quality and you want to preserve shorter reads, you can lower `length_required`.
- If the data are noisy, you can increase `qualified_quality_phred` to values like 30 or 35.
- For long-read data, it is usually safer to use a larger minimum length, for example `500`.

### Practical recommendations

- For standard Illumina data: keep `length_required` at `50` and `qualified_quality_phred` at `30`.
- For lower-quality samples: try `20` or `25` for quality and `30` or `40` for minimum length.

## 3. Choosing the contig assembler

The assembler choice is controlled in the `tool.denovo` section of [config/config.yaml](../config/config.yaml).

Example:

```yaml
tool:
  denovo:
    - 'spades'
    # - 'megahit'
    # - 'flye'
```

### When to use each option

- `spades`: best option for Illumina metagenomic and metatranscriptomic reads. It is the current default and tends to recover viral contigs well.
- `megahit`: faster and useful when you prioritize speed over maximum assembly quality.
- `flye`: recommended for long-read data (Nanopore/PacBio). Not the best choice for short Illumina reads.
- `raven`: another long-read assembler, usually simpler and faster, but with quality trade-offs.

## 4. Tuning SPAdes

If you choose SPAdes, there are two important settings to adjust:

### 4.1 Algorithm mode

In the `spades.algorithm` section you define the assembly mode:

```yaml
spades:
  algorithm:
    - '--meta'
```

- `--meta`: recommended for metagenomic data.
- `--isolate`: can work better for isolate samples or high-coverage datasets.

### 4.2 K-mer sizes

The `spades.kmer` section defines the k-mer sizes used by SPAdes:

```yaml
spades:
  kmer:
    - 'auto'
```

- `auto`: the simplest option and usually a good starting point.
- Explicit values like `21`, `33`, `55` or combinations such as `21_33_55` can be used for more specific tests.

## 5. Long-read assembly settings

If you are working with long-read data, prefer `flye` or `raven` instead of SPAdes.

In [config/config.yaml](../config/config.yaml), you can set the Flye read type:

```yaml
flye:
  type: "nano-corr"
```

You can change it to one of:

- `pacbio-raw`
- `pacbio-corr`
- `pacbio-hifi`
- `nano-raw`
- `nano-corr`
- `nano-hq`

## 6. Adjusting downstream parameters

Besides fastp and the assembler, other parameters can significantly affect results:

- `diamond.min_contig_len`: minimum contig length for DIAMOND annotation.
- `diamond.evalue`: e-value threshold for alignments.
- `kraken2.confidence`: classification confidence threshold for Kraken2.
- `palm_annot.minscore`: score threshold for RdRp candidate detection.

## 7. Recommended workflow to start

1. Start with the default values in [config/config.yaml](../config/config.yaml).
2. If the read quality is poor, tune fastp first.
3. If you have short Illumina reads, test `spades` first.
4. If you have long-read data, use `flye` or `raven`.
5. Compare results and change one parameter at a time.

## 8. Example run after editing configuration

```bash
bash NovaVir.sh --input <DIR> --output <DIR> --profile local --diamond --darkmatter
```

If you are running on a cluster, change the profile to `profile_slurm`.
