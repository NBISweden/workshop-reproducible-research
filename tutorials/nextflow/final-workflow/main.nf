#!/usr/bin/env nextflow

// This is one possible variant of the final workflow after finishing all of the
// Nextflow tutorials, not including extra material.

// Include subworkflows
include { QUALITY_CONTROLS } from "./subworkflows/quality_controls.nf"
include { ALIGNMENT        } from "./subworkflows/alignment.nf"

// Main workflow
workflow {

    main:
    // Get input files from a samplesheet
    ch_input = channel
        .fromPath ( params.input )
        .splitCsv ( header: true )

    // Define the workflow from a combination of subworkflows and processes
    DOWNLOAD_FASTQ_FILES (
        ch_input
    )
    QUALITY_CONTROLS (
        DOWNLOAD_FASTQ_FILES.out
    )
    ALIGNMENT (
        params.genome_fasta,
        DOWNLOAD_FASTQ_FILES.out
    )
    GENERATE_COUNTS_TABLE (
        ALIGNMENT.out.bam.collect(),
        params.genome_gff3
    )

    // Publish only the outputs of interest to the user
    publish:
    fastq          = DOWNLOAD_FASTQ_FILES.out
    multiqc_html   = QUALITY_CONTROLS.out.html
    multiqc_stats  = QUALITY_CONTROLS.out.general_stats
    bam            = ALIGNMENT.out.bam
    counts         = GENERATE_COUNTS_TABLE.out.counts
    counts_summary = GENERATE_COUNTS_TABLE.out.summary
}

// Output definitions
output {
    fastq {
        path "fastq"
    }
    multiqc_html {
        path "qc"
    }
    multiqc_stats {
        path "qc"
    }
    bam {
        path "bam"
    }
    counts {
        path "counts"
    }
    counts_summary {
        path "counts"
    }
}

process DOWNLOAD_FASTQ_FILES {

    // Download a single-read FASTQ file from the SciLifeLab Figshare remote

    tag "${sra_id}"
    input:
    tuple val(sra_id), val(figshare_link)

    output:
    tuple val(sra_id), path("*.fastq.gz")

    script:
    """
    curl -L -A "Mozilla/5.0" ${figshare_link} -o ${sra_id}.fastq.gz
    """
}

process GENERATE_COUNTS_TABLE {

    // Generate a count table using featureCounts.

    input:
    path(bam)
    path(annotation)

    output:
    path("counts.tsv"), emit: counts
    path("counts.tsv.summary"), emit: summary

    script:
    """
    # The transcript name is annotated as "Name" in the GFF
    featureCounts -t exon -g Name -a ${annotation} -o counts.tsv ${bam}
    """
}
