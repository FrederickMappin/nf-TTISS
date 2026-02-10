# =============================================================================
# nf-TTISS Pipeline — Multi-stage Docker Image
# =============================================================================
# This image bundles every tool the Nextflow pipeline requires so that a single
# container can be used instead of per-process biocontainer images.  It is also
# useful for local development and CI testing outside of Nextflow.
#
# Stage 1 — bioinformatics tools (FastQC, MultiQC, Cutadapt, BWA, SAMtools,
#            BEDtools) installed via Mambaforge/Conda.
# Stage 2 — slim runtime image with only the installed packages and the TTIS
#            Python application copied from stage 1.
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1: Build / install all dependencies
# ---------------------------------------------------------------------------
FROM condaforge/mambaforge:24.3.0-0 AS builder

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install bioinformatics tools into the base conda environment
RUN mamba install -y -c bioconda -c conda-forge \
        fastqc=0.12.1 \
        multiqc=1.21 \
        cutadapt=4.6 \
        bwa=0.7.17 \
        samtools=1.19 \
        bedtools=2.31.1 \
        biopython=1.83 \
        numpy=1.26 \
    && mamba clean --all --yes

# ---------------------------------------------------------------------------
# Stage 2: Slim runtime
# ---------------------------------------------------------------------------
FROM condaforge/mambaforge:24.3.0-0 AS runtime

LABEL maintainer="Frederick Mappin"
LABEL description="nf-TTISS: Targeted Transgene Integration Site Sequencing pipeline"
LABEL version="1.0.0"

# Copy the fully-resolved conda environment from the builder
COPY --from=builder /opt/conda /opt/conda

# Ensure conda binaries are on PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Create a working directory
WORKDIR /pipeline

# Copy the TTIS Python application
COPY TTIS/ /pipeline/TTIS/

# Copy the Nextflow pipeline files
COPY main.nf nextflow.config /pipeline/

# Default entrypoint — drop into bash so the image works with Nextflow
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["bash"]
