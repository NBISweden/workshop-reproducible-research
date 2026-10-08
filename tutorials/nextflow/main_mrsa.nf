#!/usr/bin/env nextflow

workflow {

    // Workflow for generating count data for the MRSA case study

    main:
    // Get input files from a samplesheet
    ch_input = channel
        .fromPath ( "samplesheet.csv" )
        .splitCsv ( header: true)

    // Define the workflow
    DOWNLOAD_FASTQ_FILES (
        ch_input
    )
    RUN_FASTQC (
        DOWNLOAD_FASTQ_FILES.out
    )
    RUN_MULTIQC (
        RUN_FASTQC.out[1].collect()
    )
    GET_GENOME_FASTA ()
    INDEX_GENOME (
        GET_GENOME_FASTA.out.fasta
    )
    ALIGN_TO_GENOME (
        DOWNLOAD_FASTQ_FILES.out,
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

    // Define which outputs to publish
    publish:
    fastq          = DOWNLOAD_FASTQ_FILES.out
    fastqc_html    = RUN_FASTQC.out[0]
    fastqc_zip     = RUN_FASTQC.out[1]
    multiqc_html   = RUN_MULTIQC.out.html
    multiqc_stats  = RUN_MULTIQC.out.general_stats
    genome_fasta   = GET_GENOME_FASTA.out.fasta
    genome_index   = INDEX_GENOME.out.index
    bam            = ALIGN_TO_GENOME.out.bam
    sorted_bam     = SORT_BAM.out.bam
    genome_gff3    = GET_GENOME_GFF3.out.gff
    counts         = GENERATE_COUNTS_TABLE.out.counts
    counts_summary = GENERATE_COUNTS_TABLE.out.summary
}

output {

    // Define where and how to publish each output

    fastq {
        path "data"
        mode "copy"
    }
    fastqc_html {
        path "fastqc"
        mode "copy"
    }
    fastqc_zip {
        path "fastqc"
        mode "copy"
    }
    multiqc_html {
        path "multiqc"
        mode "copy"
    }
    multiqc_stats {
        path "multiqc"
        mode "copy"
    }
    genome_fasta {
        path "data/ref"
        mode "copy"
    }
    genome_index {
        path "bowtie2"
        mode "copy"
    }
    bam {
        path "bam"
        mode "copy"
    }
    sorted_bam {
        path "bam"
        mode "copy"
    }
    genome_gff3 {
        path "data/ref"
        mode "copy"
    }
    counts {
        path "tables"
        mode "copy"
    }
    counts_summary {
        path "tables"
        mode "copy"
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
    curl -L ${figshare_link} -o ${sra_id}.fastq.gz
    """
}

process RUN_FASTQC {

    // Run FastQC on a FASTQ file.

    tag "${sample}"

    input:
    tuple val(sample), path(fastq)

    output:
    path("*.html")
    path("*.zip")

    script:
    """
    fastqc ${fastq} -q
    """
}

process RUN_MULTIQC {

    // Aggregate all FastQC reports into a MultiQC report.

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

    output:
    path("*.fa.gz"), emit: fasta

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O NCTC8325.fa.gz
    """
}

process GET_GENOME_GFF3 {

    // Retrieve annotation in gff3 format for a genome.

    output:
    path("*.gff3.gz"), emit: gff

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz -O NCTC8325.gff3.gz
    """
}

process INDEX_GENOME {

    // Index a genome using Bowtie 2.

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
