# DiscoVir HPC (Slurm) Usage Guide

This document explains how to execute the DiscoVir pipeline in a High-Performance Computing (HPC) environment managed by the Slurm workload manager.

## Resource Management

When running on an HPC cluster, it is crucial to use a Snakemake profile so that each pipeline rule receives the correct CPU, memory, and time limits. The DiscoVir repository provides a pre-configured profile for Slurm located at `profiles/profile_slurm`.

**The `--profile` flag is key** to running DiscoVir successfully on an HPC. It delegates the job submission to the Slurm manager, individually limiting and scheduling the partition requirements per rule.

## Basic Execution

A typical invocation for the complete deep viral discovery pipeline (Assembly + Diamond + Dark Matter) is:

```shell
bash DiscoVir.sh \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profiles/profile_slurm \
    --diamond \
    --darkmatter \
    --diamond_db /path/to/db/diamond/nr.dmnd \
    --kraken2 /path/to/db/kraken2 \
    --jobs 20
```

### Explanation of Arguments:
- `--profile profiles/profile_slurm`: Instructs Snakemake to use the Slurm executor.
- `--jobs 20`: The maximum number of parallel jobs that Snakemake will submit to the Slurm queue simultaneously.
- `--diamond` and `--darkmatter`: Enables the assembly and viral identification modules.

## Optional Modules

If you want to add the optional RVDB validation on top of the dark matter analysis, simply append the `--rvdb` flag:

```shell
bash DiscoVir.sh \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profiles/profile_slurm \
    --diamond \
    --darkmatter \
    --rvdb \
    --diamond_db /path/to/db/diamond/nr.dmnd \
    --kraken2 /path/to/db/kraken2 \
    --jobs 20
```

## Quick Profiling

If you only want a rapid initial evaluation of your samples directly on the reads without running resource-intensive assemblies, you can run:

```shell
bash DiscoVir.sh \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profiles/profile_slurm \
    --reads-kraken \
    --reads-diamond \
    --kraken2 /path/to/db/kraken2 \
    --diamond_db /path/to/db/diamond/nr.dmnd
```

For more advanced configuration of the execution parameters (e.g., K-mers, assembly algorithms), see the [Advanced Configuration Tutorial](advanced-configuration-tutorial.md).
