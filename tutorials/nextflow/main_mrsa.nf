#!/usr/bin/env nextflow

// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2

workflow {

    // Workflow for generating count data for the MRSA case study

    // Define SRA input data channel
    ch_sra_ids = Channel.fromList( ["SRR935090", "SRR935091", "SRR935092"] )

    // Define the workflow
    GET_SRA_BY_ACCESSION (
        ch_sra_ids
    )
    RUN_FASTQC (
        GET_SRA_BY_ACCESSION.out
    )
    RUN_MULTIQC (
        RUN_FASTQC.out[1].collect()
    )
    GET_GENOME_FASTA ()
    INDEX_GENOME (
        GET_GENOME_FASTA.out.fasta
    )
    ALIGN_TO_GENOME (
        GET_SRA_BY_ACCESSION.out,
        INDEX_GENOME.out.index
    )
    SORT_BAM (
        ALIGN_TO_GENOME.out.bam
    )
    GET_GENOME_GFF3 ()
    GENERATE_COUNTS_TABLE (
        SORT_BAM.out.bam.collect(),
        GET_GENOME_GFF3.out.gff
    )
}

process GET_SRA_BY_ACCESSION {

    // Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run
    // accession number.

    tag "${sra_id}"
    publishDir "results/data/raw_internal",
        mode: "copy"

    input:
    val(sra_id)

    output:
    tuple val(sra_id), path("*.fastq.gz")

    script:
    """
    fastq-dump ${sra_id} \
        -X 25000 \
        --readids \
        --dumpbase \
        --skip-technical \
        --gzip \
        -Z \
        > ${sra_id}.fastq.gz
    """
}

process RUN_FASTQC {

    // Run FastQC on a FASTQ file.

    tag "${sample}"
    publishDir "results/", mode: "copy"

    input:
    tuple val(sample), path(fastq)

    output:
    path("*.html")
    path("*.zip")

    script:
    """
    # Run FastQC
    fastqc ${fastq} -q
    """
}

process RUN_MULTIQC {

    // Aggregate all FastQC reports into a MultiQC report.

    publishDir "results/results",
        pattern: "*.html",
        mode: "copy"
    publishDir "results/intermediate",
        pattern: "*.txt",
        mode: "copy"

    input:
    path(zips)

    output:
    path("*.html"), emit: html
    path("multiqc_general_stats.txt"), emit: general_stats

    script:
    """
    multiqc -n multiqc.html .
    mv multiqc_data/multiqc_general_stats.txt multiqc_general_stats.txt
    """
}

process GET_GENOME_FASTA {

    // Retrieve the sequence in fasta format for a genome.

    publishDir "results/data/raw_external",
        mode: "copy"

    output:
    path("*.fa.gz"), emit: fasta

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O NCTC8325.fa.gz
    """
}

process GET_GENOME_GFF3 {

    // Retrieve annotation in gff3 format for a genome.

    publishDir "results/data/raw_external",
        mode: "copy"

    output:
    path("*.gff3.gz"), emit: gff

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz -O NCTC8325.gff3.gz
    """
}

process INDEX_GENOME {

    // Index a genome using Bowtie 2.

    publishDir "results/intermediate/",
        mode: "copy"

    input:
    path(fasta)

    output:
    path("*.bt2"), emit: index

    script:
    """
    # Bowtie2 cannot use .gz, so unzip to a temporary file first
    gunzip -c ${fasta} > tempfile
    bowtie2-build tempfile NCTC8325
    """
}

process ALIGN_TO_GENOME {

    // Align a fastq file to a genome index using Bowtie 2.

    tag "${sample}"
    publishDir "results/intermediate/",
        mode: "copy"

    input:
    tuple val(sample), path(fastq)
    path(idx)

    output:
    tuple val(sample), path("*.bam"), emit: bam

    script:
    """
    bowtie2 -x NCTC8325 -U ${fastq} > ${sample}.bam
    """
}

process SORT_BAM {

    // Sort a bam file.

    tag "${sample}"
    publishDir "results/intermediate/",
        mode: "copy"

    input:
    tuple val(sample), path(bam)

    output:
    path("*.sorted.bam"), emit: bam

    script:
    """
    samtools sort ${bam} > ${sample}.sorted.bam
    """
}

process GENERATE_COUNTS_TABLE {

    // Generate a count table using featureCounts.

    publishDir "results/results/tables",
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
