# TTISS Pipeline - Targeted Transgene Integration Site Sequencing

A Nextflow pipeline for analyzing off-target integration sites in CRISPR-based genome editing experiments using dsODN donors.

## Overview

This pipeline processes paired-end sequencing reads to identify and characterize genome integration sites, performing quality control, adapter trimming, read mapping, window extraction, and guide RNA matching.

## Pipeline Steps

1. **FastQC** - Quality control of raw sequencing reads
2. **MultiQC** - Aggregate FastQC reports into a single summary
3. **CUTADAPT** - Filter reads by UMI primer pattern (`^NNNNNNNNNNNNNNNNNCGC`)
4. **TRIM_AND_TRUNCATE** - Trim primers and truncate reads (R1: 25bp, R2: 15bp)
5. **BWA_MAP** - Map reads to reference genome with ultra-sensitive parameters
6. **EXTRACT_WINDOWS** - Extract 100bp windows around integration sites (≥2 reads after merging)
7. **MATCH_GUIDES** - Match guide RNA sequences to extracted windows

## Requirements

- [Nextflow](https://www.nextflow.io/) (≥23.0.0)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/FrederickMappin/nf-TTISS.git
cd nf-TTISS
```

### 2. Prepare your data

**Required files:**
- Paired-end FASTQ files (gzipped): `data/*_R{1,2}.fastq.gz`
- Reference genome (FASTA): `data/reference/genome.fa`
- BWA index files for reference (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`)
- Guide RNA file (CSV): `TTIS/Guides.txt`

**Guide file format:**
```
1,GAGTCCGTGAGAGGCAGAGG
2,AATTCTTCTTCTGAAATAGT
```

### 3. Run the pipeline

**Basic usage:**
```bash
nextflow run main.nf -profile docker
```

**With custom parameters:**
```bash
nextflow run main.nf -profile docker \
  --reads "data/*_R{1,2}.fastq.gz" \
  --reference "data/reference/genome.fa" \
  --guides "TTIS/Guides.txt" \
  --pam "NGG" \
  --outdir "results"
```

**Resume a previous run:**
```bash
nextflow run main.nf -profile docker --resume
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | `data/*_R{1,2}.fastq.gz` | Input paired-end FASTQ files |
| `--outdir` | `results` | Output directory |
| `--reference` | `data/reference/hg19_chr8.fa` | Reference genome (FASTA) |
| `--guides` | `TTIS/Guides.txt` | Guide RNA sequences (CSV) |
| `--pam` | `NGG` | PAM sequence for guide matching |
| `--primer` | `^NNNNNNNNNNNNNNNNNCGC` | UMI primer pattern (17bp + CGC) |
| `--error_rate` | `0.1` | Cutadapt error rate |
| `--min_length` | `40` | Minimum read length after trimming |
| `--quality_cutoff` | `20` | Quality score cutoff |

## Output Structure

```
results/
├── fastqc/              # FastQC quality control reports
├── multiqc/             # MultiQC aggregate report
│   └── multiqc_report.html
├── filtered/            # Primer-filtered reads
├── trimmed/             # Trimmed and truncated reads
├── mapped/              # BAM files from BWA mapping
├── windows/             # Extracted 100bp windows
│   ├── *_windows.bed    # Window coordinates
│   └── *_windows.fasta  # Window sequences
├── matches/             # Guide RNA matching results
│   └── *_matches.csv
└── reports/             # Pipeline execution reports
    ├── execution_report.html
    ├── timeline.html
    └── trace.txt
```

## Key Features

- **UMI-based filtering**: Retains reads with correct primer pattern
- **Ultra-sensitive mapping**: BWA parameters optimized for short reads (-k 10 -T 15)
- **Window merging**: Combines integration sites within 10bp
- **Strand-aware matching**: Searches both forward and reverse complement
- **Reproducible**: Docker/Singularity containers ensure consistency

