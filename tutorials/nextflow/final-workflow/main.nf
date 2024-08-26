#!/use/bin/env Nextflow

// This is one possible variant of the final workflow after finishing all of the
// Nextflow tutorials, not including extra material.

// Include subworkflows
include { QUALITY_CONTROLS } from "./subworkflows/quality_controls.nf"
include { ALIGNMENT        } from "./subworkflows/alignment.nf"

// Main workflow
workflow {

    // Get input files from a samplesheet
    ch_input = Channel
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
        ALIGNMENT.out.bam,
        params.genome_gff3
    )
}

process DOWNLOAD_FASTQ_FILES {

    // Download a single-read FASTQ file from the SciLifeLab Figshare remote

    tag "${sra_id}"
    publishDir "${params.outdir}/fastq",
        mode: "copy"

    input:
    tuple val(sra_id), val(figshare_link)

    output:
    tuple val(sra_id), path("*.fastq.gz")

    script:
    """
    wget ${figshare_link} -O ${sra_id}.fastq.gz
    """
}

process GENERATE_COUNTS_TABLE {

    // Generate a count table using featureCounts.

    publishDir "${params.outdir}/counts",
        mode: "copy"

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
