#!/usr/bin/env nextflow

// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2

// Define where results should end up
resultsdir = "results/"

workflow {

    // Workflow for generating count data for the MRSA case study

    // Define SRA input data channel
    Channel
        .fromList( params.sra_id_list )
        .set{ sra_ids }

    // Define the workflow
    get_sra_by_accession(sra_ids)
    run_fastqc(get_sra_by_accession.out.sra_data)
    run_multiqc(run_fastqc.out.zip.collect())
    get_genome_fasta()
    index_genome(get_genome_fasta.out.fasta)
    align_to_genome(get_sra_by_accession.out.sra_data,
                    index_genome.out.index)
    sort_bam(align_to_genome.out.bam)
    get_genome_gff3()
    generate_counts_table(sort_bam.out.bam.collect(),
                          get_genome_gff3.out.gff)
}

process get_sra_by_accession {

    // Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run
    // accession number.

    tag "${sra_id}"
    publishDir "${resultsdir}/data/",
        mode: "copy"

    input:
    val(sra_id)

    output:
    tuple val(sra_id), path("${sra_id}.fastq.gz"), emit: sra_data

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

process run_fastqc {

    // Run FastQC on a FASTQ file.

    tag "${sample}"
    publishDir "${resultsdir}/qc",
        mode: "copy",
        pattern: "*.html"
    publishDir "${resultsdir}/qc/intermediate",
        mode: "copy",
        pattern: "*.zip"

    input:
    tuple val(sample), path(fastq)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:
    """
    # Run FastQC
    fastqc ${fastq} -q
    """
}

process run_multiqc {

    // Aggregate all FastQC reports into a MultiQC report.

    publishDir "${resultsdir}/qc/",
        mode: "copy"

    input:
    path(zips)

    output:
    path("*.html"), emit: html
    path("multiqc_data/multiqc_general_stats.txt"), emit: general_stats

    script:
    """
    multiqc -n multiqc.html .
    """
}

process get_genome_fasta {

    // Retrieve the sequence in fasta format for a genome.

    publishDir "${resultsdir}/idx/",
        mode: "copy"

    output:
    path("*.fa.gz"), emit: fasta

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
    """
}

process get_genome_gff3 {

    // Retrieve annotation in gff3 format for a genome.

    publishDir "${resultsdir}/idx/",
        mode: "copy"

    output:
    path("*.gff3.gz"), emit: gff

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz
    """
}

process index_genome {

    // Index a genome using Bowtie 2.

    publishDir "${resultsdir}/idx/",
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

process align_to_genome {

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

process sort_bam {

    // Sort a bam file.

    tag "${sample}"
    publishDir "${resultsdir}/bam/",
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

process generate_counts_table {

    // Generate a count table using featureCounts.

    publishDir "${resultsdir}/",
        mode: "copy"

    input:
    path(bam)
    path(annotation)

    output:
    path("counts.tsv"), emit: counts

    script:
    """
    # The transcript name is annotated as "Name" in the GFF
    featureCounts -t exon -g Name -a ${annotation} -o counts.tsv ${bam}
    """
}
