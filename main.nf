#!/usr/bin/env nextflow

/*
 * Nextflow pipeline with Cutadapt
 * Adapter trimming workflow
 */

nextflow.enable.dsl=2

// Parameters
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"
params.primer = "^NNNNNNNNNNNNNNNNNCGC"
params.error_rate = 0.1
params.min_length = 40
params.quality_cutoff = 20
params.reference = "data/reference/hg19_chr8.fa"
params.guides = "TTIS/Guides.txt"
params.pam = "NGG"

// Print parameters
log.info """\
    CUTADAPT PIPELINE - PRIMER FILTERING
    =====================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    primer       : ${params.primer}
    error_rate   : ${params.error_rate}
    min_length   : ${params.min_length}
    quality      : ${params.quality_cutoff}
    """
    .stripIndent()

/*
 * Process: Run FastQC for quality control
 */
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.{html,zip}", emit: reports
    
    script:
    """
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

/*
 * Process: Run MultiQC to aggregate FastQC reports
 */
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path('*')
    
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    
    script:
    """
    multiqc .
    """
}

/*
 * Process: Run Cutadapt for primer filtering
 */
process CUTADAPT {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered_R*.fastq.gz"), emit: filtered_reads
    path "${sample_id}_cutadapt_report.txt", emit: log
    
    script:
    """
    cutadapt \
        -g "${params.primer}" \
        -e ${params.error_rate} \
        --discard-untrimmed \
        --minimum-length ${params.min_length} \
        --quality-cutoff ${params.quality_cutoff} \
        --rename='{id} umi:{cut_prefix}' \
        --action=none \
        -o ${sample_id}_filtered_R1.fastq.gz \
        -p ${sample_id}_filtered_R2.fastq.gz \
        ${reads[0]} ${reads[1]} \
        > ${sample_id}_cutadapt_report.txt
    """
}

/*
 * Process: Trim primer and truncate to mapping lengths
 */
process TRIM_AND_TRUNCATE {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R*.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_trim_report.txt", emit: log
    
    script:
    """
    # Step 1: Trim the primer from R1 and R2
    cutadapt \
        -g "${params.primer}" \
        -o ${sample_id}_primer_trimmed_R1.fastq.gz \
        -p ${sample_id}_primer_trimmed_R2.fastq.gz \
        ${reads[0]} ${reads[1]} \
        > ${sample_id}_trim_report.txt
    
    # Step 2: Truncate R1 to 25bp
    cutadapt \
        --length 25 \
        -o ${sample_id}_trimmed_R1.fastq.gz \
        ${sample_id}_primer_trimmed_R1.fastq.gz \
        >> ${sample_id}_trim_report.txt
    
    # Step 3: Truncate R2 to 15bp
    cutadapt \
        --length 15 \
        -o ${sample_id}_trimmed_R2.fastq.gz \
        ${sample_id}_primer_trimmed_R2.fastq.gz \
        >> ${sample_id}_trim_report.txt
    """
}

/*
 * Process: Map trimmed reads with BWA
 */
process BWA_MAP {
    tag "$sample_id"
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference
    path reference_index
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.bam.bai"), emit: bai
    
    script:
    """
    bwa mem -t ${task.cpus} -k 10 -T 15 -A 1 -B 1 -O 1 -E 1 -I 1000 ${reference} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    
    samtools view -f 2 -Sb ${sample_id}.sam | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    
    samtools index ${sample_id}.bam
    """
}

/*
 * Process: Extract 100bp windows around integration sites
 */
process EXTRACT_WINDOWS {
    tag "$sample_id"
    publishDir "${params.outdir}/windows", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}_windows.bed"), emit: bed
    tuple val(sample_id), path("${sample_id}_windows.fasta"), emit: fasta
    
    script:
    """
    # Install tools if needed
    mamba install -y -c bioconda -c conda-forge samtools bedtools
    
    # Create FASTA index if it doesn't exist
    samtools faidx ${reference}
    
    # Create windows with read counts (no filtering yet)
    samtools view -F 4 ${bam} | \
      awk 'BEGIN{OFS="\\t"} {print \$3, \$4-50, \$4+50}' | \
      sort -k1,1 -k2,2n | \
      uniq -c | \
      awk '{print \$2"\\t"\$3"\\t"\$4"\\t"\$1}' > ${sample_id}_windows_raw.bed
    
    # Merge windows within 10bp and sum read counts, then filter for ≥2 total reads
    bedtools merge -i ${sample_id}_windows_raw.bed -c 4 -o sum -d 10 | \
      awk '\$4 >= 1' > ${sample_id}_windows.bed
    
    # Extract sequences from reference
    bedtools getfasta -fi ${reference} -bed ${sample_id}_windows.bed \
      -fo ${sample_id}_windows.fasta
    """
}

/*
 * Process: Match guides to extracted windows
 */
process MATCH_GUIDES {
    tag "$sample_id"
    publishDir "${params.outdir}/matches", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fasta)
    path guides
    
    output:
    tuple val(sample_id), path("${sample_id}_matches.csv"), emit: matches
    
    script:
    """
    # Install Python dependencies
    pip3 install biopython numpy
    
    # Run the matcher
    python3 ${projectDir}/TTIS/main.py \
        --guidefile ${guides} \
        --fastafile ${fasta} \
        --pam ${params.pam} \
        --output ${sample_id}_matches.csv
    """
}

/*
 * Main workflow
 */
workflow {
    // Create channel from input files
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
    
    // Reference genome and index files
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    reference_index_ch = Channel.fromPath("${params.reference}.*", checkIfExists: true).collect()
    
    // Guide file
    guides_ch = Channel.fromPath(params.guides, checkIfExists: true)
    
    // Run FastQC on raw reads
    FASTQC(read_pairs_ch)
    
    // Run MultiQC to aggregate FastQC reports
    MULTIQC(FASTQC.out.reports.collect())
    
    // Run Cutadapt to filter by primer
    CUTADAPT(read_pairs_ch)
    
    // Trim primer and truncate to mapping lengths
    TRIM_AND_TRUNCATE(CUTADAPT.out.filtered_reads)
    
    // Map trimmed reads with BWA
    BWA_MAP(TRIM_AND_TRUNCATE.out.trimmed_reads, reference_ch, reference_index_ch)
    
    // Extract 100bp windows around integration sites with ≥2 reads
    EXTRACT_WINDOWS(BWA_MAP.out.bam, reference_ch)
    
    // Match guides to windows
    MATCH_GUIDES(EXTRACT_WINDOWS.out.fasta, guides_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Results are in ${params.outdir}\n" : "Oops .. something went wrong" )
}
