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
 * Main workflow
 */
workflow {
    // Create channel from input files
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
    
    // Run Cutadapt to filter by primer
    CUTADAPT(read_pairs_ch)
    
    // Trim primer and truncate to mapping lengths
    TRIM_AND_TRUNCATE(CUTADAPT.out.filtered_reads)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Results are in ${params.outdir}\n" : "Oops .. something went wrong" )
}
