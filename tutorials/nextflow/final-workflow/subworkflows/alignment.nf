//
// Alignment subworkflow
//
workflow ALIGNMENT {
    
    take:
    fasta
    fastq

    main:
    INDEX_GENOME (
        fasta
    )
    ALIGN_TO_GENOME (
        fastq,
        INDEX_GENOME.out.index
    )
    SORT_BAM (
        ALIGN_TO_GENOME.out.bam
    )

    emit:
    bam = SORT_BAM.out.bam
}

process INDEX_GENOME {

    // Index a genome using Bowtie 2.
    //
    // Publishing is not needed, as output only contain index files used by
    // other processes and do not need to be inspected by the user.

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
    //
    //  Publishing is not needed, as the subsequent process produces a sorted
    //  BAM that the user can view if desired.

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
    publishDir "${params.outdir}/bam/",
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
