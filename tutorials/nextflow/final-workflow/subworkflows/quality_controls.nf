//
// Quality controls subworkflow
//
workflow QUALITY_CONTROLS {

    take:
    fastq

    main:
    RUN_FASTQC (
        fastq
    )
    RUN_MULTIQC (
        RUN_FASTQC.out.zip.collect()
    )

    emit:
    html          = RUN_MULTIQC.out.html
    general_stats = RUN_MULTIQC.out.general_stats
}

process RUN_FASTQC {

    // Run FastQC on a FASTQ file.
    //
    // No output publishing required as all results are aggregated by MultiQC

    tag "${sample}"

    input:
    tuple val(sample), path(fastq)

    output:
    path("*.html"), emit: html
    path("*.zip"),  emit: zip

    script:
    """
    fastqc ${fastq} -q
    """
}

process RUN_MULTIQC {

    // Aggregate all FastQC reports into a MultiQC report.

    publishDir "${params.outdir}/qc",
        mode: "copy"

    input:
    path(zips)

    output:
    path("*.html")                   , emit: html
    path("multiqc_general_stats.txt"), emit: general_stats

    script:
    """
    multiqc -n multiqc.html .
    mv multiqc_data/multiqc_general_stats.txt multiqc_general_stats.txt
    """
}
