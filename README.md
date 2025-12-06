# TTISS - Nextflow Pipeline with Cutadapt

A Nextflow pipeline for filtering reads by primer sequence using Cutadapt. Only keeps reads that start with a specific primer pattern and extracts UMI sequences.

## Features

- **Primer filtering**: Keep only reads starting with specified primer
- **UMI extraction**: Extract and add UMI to read names from primer sequence
- **Quality filtering**: Remove low-quality bases
- **Discard untrimmed**: Only keep reads with the primer
- Multiple execution profiles (local, Docker, Singularity, SLURM)

## Requirements

- Nextflow >= 23.0.0
- Java 11 or later
- Docker, Singularity, or Conda (depending on profile)

## Installation

Nextflow is already installed. To verify:

```bash
nextflow -version
```

## Quick Start

1. **Prepare your data**: Place paired-end FASTQ files in a `data/` directory
   - Files should follow the pattern: `*_R1.fastq.gz` and `*_R2.fastq.gz`

2. **Run the pipeline**:

```bash
# Local execution (requires tools installed)
nextflow run main.nf

# Using Docker
nextflow run main.nf -profile docker

# Using Singularity
nextflow run main.nf -profile singularity

# On SLURM cluster
nextflow run main.nf -profile slurm
```

## Configuration

### Input Parameters

Edit `nextflow.config` or provide parameters via command line:

```bash
nextflow run main.nf \
    --reads "data/*_R{1,2}.fastq.gz" \
    --outdir "results" \
    --primer "^NNNNNNNNNNNNNNNNNCGC" \
    --error_rate 0.1 \
    --min_length 40 \
    --quality_cutoff 20
```

### Parameters Explained

- `--primer`: Primer pattern to search for at the start of reads (^ = anchored at 5' end)
  - `N` = any nucleotide (UMI positions)
  - Use specific bases for the fixed primer sequence
- `--error_rate`: Maximum allowed error rate (0.1 = 10%)
- `--min_length`: Minimum read length after trimming
- `--quality_cutoff`: Quality threshold for 3' end trimming

### Profiles

- **standard**: Local execution
- **docker**: Run with Docker containers
- **singularity**: Run with Singularity containers
- **slurm**: Submit jobs to SLURM scheduler

## Output

Results are saved to `results/` directory:

```
results/
├── filtered/             # Filtered FASTQ files and Cutadapt reports
└── reports/              # Pipeline execution reports
    ├── execution_report.html
    ├── timeline.html
    ├── trace.txt
    └── dag.html
```

### Understanding the Output

- **Filtered reads**: Only reads starting with the primer are kept
- **UMI in read names**: The primer sequence is added to read names as `umi:SEQUENCE`
- **Cutadapt report**: Shows how many reads passed filtering

## Example

```bash
# Create test data directory
mkdir -p data

# Run with custom parameters
nextflow run main.nf \
    --reads "data/sample*_R{1,2}.fastq.gz" \
    --outdir "my_results" \
    --quality_cutoff 30 \
    -profile docker
```

## Resume

If the pipeline fails, you can resume from the last successful step:

```bash
nextflow run main.nf -resume
```

## Clean Up

Remove work directory and logs:

```bash
rm -rf work/ .nextflow/ .nextflow.log*
```

## Support

For issues or questions, please refer to the Nextflow documentation:
- https://www.nextflow.io/docs/latest/
- https://cutadapt.readthedocs.io/