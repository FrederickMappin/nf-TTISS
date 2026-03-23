# nf-TTISS —  tagmentation-based tag integration site

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A523.0.0-brightgreen?logo=nextflow)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-ready-blue?logo=docker)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

[Highly Parallel Profiling of Cas9 Variant Specificity]
A  Nextflow pipeline Tagmentation-based Tag Integration Site Sequencing (TTISS), a
pipeline for analyzing double-strand breaks, such as those created by CRISPR
nucleases.
---

## Table of Contents

- [Background](#background)
- [Pipeline Overview](#pipeline-overview)
- [Architecture](#architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Docker Usage](#docker-usage)
- [Parameters Reference](#parameters-reference)
- [Input Files](#input-files)
- [Output Structure](#output-structure)
- [Guide Matching Algorithm](#guide-matching-algorithm)
- [Execution Profiles](#execution-profiles)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

---

## Background

After CRISPR-mediated genome editing with dsODN donors, the donor DNA can integrate at both the on-target cut site and at off-target sites throughout the genome.  **TTISS** (Targeted Transgene Integration Site Sequencing) identifies those sites by:

1. Sequencing fragmented genomic DNA that contains the dsODN tag.
2. Mapping the reads to a reference genome.
3. Extracting genomic windows around mapped positions.
4. Matching guide-RNA sequences to those windows to determine which CRISPR target site is responsible.

---

## Pipeline Overview

The pipeline is implemented in seven Nextflow processes that execute sequentially:

| # | Process | Tool(s) | Purpose |
|---|---------|---------|---------|
| 1 | **FASTQC** | FastQC 0.12 | Per-read quality control |
| 2 | **MULTIQC** | MultiQC 1.21 | Aggregate QC into a single HTML report |
| 3 | **CUTADAPT** | Cutadapt 4.6 | Filter reads by UMI primer pattern `^NNNNNNNNNNNNNNNNNCGC` (17 random bases + CGC anchor) |
| 4 | **TRIM_AND_TRUNCATE** | Cutadapt 4.6 | Remove primer, then truncate R1→25 bp, R2→15 bp |
| 5 | **BWA_MAP** | BWA 0.7 + SAMtools | Map reads with ultra-sensitive parameters (`-k 10 -T 15`) |
| 6 | **EXTRACT_WINDOWS** | SAMtools + BEDtools | Extract 100 bp windows around mapped positions; merge within 10 bp; keep windows with ≥1 supporting read |
| 7 | **MATCH_GUIDES** | Python (BioPython, NumPy) | Score each window against every guide RNA (+ reverse complement) using Hamming distance and PAM validation |

---

## Architecture

```
data/*_R{1,2}.fastq.gz
        │
        ▼
   ┌──────────┐    ┌──────────┐
   │  FASTQC  │───▶│ MULTIQC  │──▶ multiqc_report.html
   └──────────┘    └──────────┘
        │
        ▼
   ┌──────────┐
   │ CUTADAPT │  filter by UMI primer
   └──────────┘
        │
        ▼
  ┌──────────────────┐
  │ TRIM_AND_TRUNCATE│  trim primer + truncate
  └──────────────────┘
        │
        ▼
   ┌──────────┐
   │ BWA_MAP  │  align to reference genome
   └──────────┘
        │
        ▼
  ┌─────────────────┐
  │ EXTRACT_WINDOWS │  100 bp windows (.bed + .fasta)
  └─────────────────┘
        │
        ▼
  ┌──────────────┐
  │ MATCH_GUIDES │  guide RNA matching (.csv)
  └──────────────┘
```

---

## Requirements

| Dependency | Version | Notes |
|------------|---------|-------|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.0.0 | Workflow engine |
| [Docker](https://www.docker.com/) **or** [Singularity](https://sylabs.io/singularity/) | any recent | Container runtime |
| Java | 11–17 | Required by Nextflow |

> **Tip:** If you use the bundled **Dockerfile** all bioinformatics tools are pre-installed — you only need Nextflow + Docker on the host.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/FrederickMappin/nf-TTISS.git
cd nf-TTISS

# 2. (Optional) Build the all-in-one Docker image
docker build -t nf-ttiss:latest .

# 3. Verify Nextflow is available
nextflow -version
```

---

## Quick Start

### Prepare Your Data

Place the following files in the project directory:

| File | Location | Description |
|------|----------|-------------|
| Paired-end FASTQ files (gzipped) | `data/*_R{1,2}.fastq.gz` | Raw sequencing reads |
| Reference genome (FASTA) | `data/reference/genome.fa` | e.g. hg19, hg38, or custom |
| BWA index files | `data/reference/genome.fa.{amb,ann,bwt,pac,sa}` | Pre-built with `bwa index` |
| Guide RNA file | `TTIS/Guides.txt` | CSV: `<id>,<20-nt guide sequence>` |

**Guide file format (`TTIS/Guides.txt`):**
```
1,AATTCTTCTTCTGAAATAGT
2,AATTCTACATCTGAAATAGT
```

### Run the Pipeline

```bash
# Basic — uses per-process biocontainer images
nextflow run main.nf -profile docker

# Custom parameters
nextflow run main.nf -profile docker \
  --reads  "data/*_R{1,2}.fastq.gz" \
  --reference "data/reference/hg38.fa" \
  --guides "TTIS/Guides.txt" \
  --pam    "NGG" \
  --outdir "results"

# Resume a failed/interrupted run
nextflow run main.nf -profile docker -resume
```

---

## Docker Usage

### Option A — Per-process Containers (Default)

When you run with `-profile docker`, Nextflow automatically pulls the correct biocontainer for every process.  No manual build step required.

### Option B — Single All-in-one Image

Build the provided `Dockerfile` which bundles **every** tool (FastQC, MultiQC, Cutadapt, BWA, SAMtools, BEDtools, BioPython, NumPy):

```bash
# Build
docker build -t nf-ttiss:latest .

# Interactive shell inside the image
docker run --rm -it \
  -v "$(pwd)/data:/pipeline/data:ro" \
  -v "$(pwd)/results:/pipeline/results" \
  nf-ttiss:latest bash

# Run the full pipeline inside the container
docker run --rm \
  -v "$(pwd)/data:/pipeline/data:ro" \
  -v "$(pwd)/results:/pipeline/results" \
  nf-ttiss:latest \
  "nextflow run /pipeline/main.nf -profile standard"
```

### Option C — Docker Compose

```bash
# Build the image
docker compose build

# Drop into a shell
docker compose run --rm pipeline bash

# Run only the Python guide-matcher
docker compose run --rm ttis-matcher \
  --guidefile TTIS/Guides.txt \
  --fastafile data/windows.fasta \
  --pam NGG
```

---

## Parameters Reference

### Input / Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `data/*_R{1,2}.fastq.gz` | Glob pattern for paired-end FASTQ files |
| `--outdir` | `results` | Directory for all pipeline outputs |
| `--reference` | `data/reference/hg19_chr8.fa` | Reference genome in FASTA format |
| `--guides` | `TTIS/Guides.txt` | Guide RNA file (CSV: `id,sequence`) |
| `--pam` | `NGG` | PAM sequence for guide matching (`N` = any nucleotide) |

### Read Processing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--primer` | `^NNNNNNNNNNNNNNNNNCGC` | UMI primer pattern — 17 random bases (UMI) followed by CGC anchor |
| `--error_rate` | `0.1` | Maximum error rate for Cutadapt primer matching |
| `--min_length` | `40` | Discard reads shorter than this after trimming |
| `--quality_cutoff` | `20` | Phred quality score cutoff for 3′ trimming |

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_memory` | `8 GB` | Upper memory bound across all processes |
| `max_cpus` | `4` | Upper CPU bound across all processes |
| `max_time` | `24 h` | Upper wall-time bound across all processes |

---

## Input Files

### Paired-end FASTQ Files

Gzipped paired-end reads following the naming convention `<sample>_R{1,2}.fastq.gz`.  Nextflow automatically groups R1/R2 pairs via `fromFilePairs`.

### Reference Genome

A FASTA file plus its BWA index.  Generate the index once with:

```bash
bwa index data/reference/genome.fa
```

This creates `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` files alongside the FASTA.

### Guide RNA File

A plain-text CSV with one guide per line:

```
<numeric_id>,<20-nucleotide guide sequence>
```

Example:
```
1,AATTCTTCTTCTGAAATAGT
2,AATTCTACATCTGAAATAGT
```

The pipeline automatically generates reverse complements of every guide for matching.

---

## Output Structure

```
results/
├── fastqc/                     # Per-sample FastQC HTML + ZIP reports
├── multiqc/
│   ├── multiqc_report.html     # Aggregated quality report
│   └── multiqc_data/           # MultiQC raw data
├── filtered/                   # Primer-filtered reads (Cutadapt)
│   ├── <sample>_filtered_R1.fastq.gz
│   └── <sample>_filtered_R2.fastq.gz
├── trimmed/                    # Trimmed + truncated reads
│   ├── <sample>_trimmed_R1.fastq.gz
│   └── <sample>_trimmed_R2.fastq.gz
├── mapped/                     # Aligned reads
│   ├── <sample>.sam
│   ├── <sample>.bam
│   └── <sample>.bam.bai
├── windows/                    # Genomic windows around integration sites
│   ├── <sample>_windows.bed    # Coordinates + read counts
│   └── <sample>_windows.fasta  # Extracted sequences
├── matches/                    # Guide RNA matching results
│   └── <sample>_matches.csv    # Columns: Name, SingleMatch, Mismatches,
│                               #   GuideNumber, CrudeSpacer, RealSpacer,
│                               #   CutSiteScore, PAM_Match
└── reports/                    # Nextflow execution reports
    ├── execution_report.html
    ├── timeline.html
    ├── trace.txt
    └── dag.html
```

### Understanding `_matches.csv`

| Column | Description |
|--------|-------------|
| `Name` | FASTA header of the genomic window |
| `SingleMatch` | `True` if only one guide matched below the mismatch threshold |
| `Mismatches` | Hamming distance between the best-matching spacer and the guide |
| `GuideNumber` | Numeric ID of the matched guide from the input file |
| `CrudeSpacer` | Raw aligned sequence with PAM context |
| `RealSpacer` | Annotated spacer: matching positions shown as `-`, mismatches as the actual nucleotide |
| `CutSiteScore` | Distance (bp) between the predicted cut site and the integration site |
| `PAM_Match` | `True` if the PAM adjacent to the spacer matches the `--pam` pattern |

---

## Guide Matching Algorithm

For each 100 bp window the matcher:

1. Slides a 20 bp window across the sequence.
2. Computes the **Hamming distance** to every guide (forward + reverse complement).
3. Keeps the guide with the lowest mismatch count (threshold ≤ 6).
4. Determines orientation (forward vs. reverse complement) and extracts the PAM-adjacent bases.
5. Validates the PAM against the user-supplied pattern (`N` = wildcard).
6. Annotates mismatches: matching positions → `-`; mismatches → actual base.

Only rows where the PAM matches are written to the output CSV.

---

## Execution Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution (no containers) — requires all tools installed locally |
| `docker` | Per-process Docker containers from BioContainers registry |
| `singularity` | Per-process Singularity images (for HPC clusters) |
| `slurm` | SLURM job scheduler (combine with `docker` or `singularity`) |

Combine profiles:
```bash
nextflow run main.nf -profile singularity,slurm
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `No such file: data/reference/genome.fa.*` | Run `bwa index data/reference/genome.fa` to create BWA index files |
| Nextflow `Process requirement exceeds available memory` | Lower `max_memory` in `nextflow.config` or increase Docker memory |
| `EXTRACT_WINDOWS` fails with `mamba: command not found` | Use `-profile docker` so the process runs inside the mambaforge container |
| Empty `_matches.csv` output | No windows matched any guide with ≤ 6 mismatches **and** a valid PAM; check your guide file and PAM setting |
| `docker: permission denied` | Add your user to the `docker` group: `sudo usermod -aG docker $USER` |

---

## Contributing

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/<name>`.
3. Commit your changes with clear messages.
4. Open a pull request against `main`.

Please ensure any new Nextflow processes have corresponding container definitions in both the `docker` and `singularity` profiles.

---

## License

This project is released under the [MIT License](LICENSE).

## Key Features

- **UMI-based filtering**: Retains reads with correct primer pattern
- **Ultra-sensitive mapping**: BWA parameters optimized for short reads (-k 10 -T 15)
- **Window merging**: Combines integration sites within 10bp
- **Strand-aware matching**: Searches both forward and reverse complement
- **Reproducible**: Docker/Singularity containers ensure consistency

