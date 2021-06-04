#!/usr/bin/env nextflow

// Enable DSL2 for a more modular workflow
nextflow.enable.dsl=2

// Define where results should end up
resultsdir = "results/"

workflow {
    // Get the input SRA IDs from parameters as channel
    Channel
        .from( params.sra_ids )
        .set{ sra_ids }

    // Define the workflow
    get_sra_by_accession(sra_ids)
    run_fastqc(get_sra_by_accession.out)
    run_multiqc(run_fastqc.out.zip.collect())
}

process get_sra_by_accession {
    tag "${sra_id}"
    publishDir "${resultsdir}/data/",
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

process run_fastqc {
    tag "${sample}"
    publishDir "${resultsdir}/",
        mode: "copy",
        // Here we specify that anything ending in `.zip` should be put into the
        // `intermediate/` directory rathern in the default publishDir. One
        // could argue, however, that since the compressed output of FastQC are
        // just intermediate files for MultiQC we do not need to publish them.
        saveAs: { filename ->
            filename.indexOf(".zip") > 0 ? \
                "${filename}" : "intermediate/${filename}"
        }

    input:
    // Note that we name the input value "sample" here, rather than sticking
    // with the original "sra_id" name given in the previous process. This works
    // since each process functions like a function, meaning that we can input
    // any parameter we want (that is named whatever we want), but it is
    // internally named something different inside the process.
    tuple val(sample), path(fastq)

    output:
    // Since the compressed files will be used by MultiQC (but not the HTML) we
    // "emit" these files and name them, which can then be used by other
    // processes specifically.
    path("*.html")
    path("*.zip"), emit: zip

    script:
    """
    # Run FastQC
    fastqc ${fastq} -q

    # Note that We don't need to specify `-o .`, as all Nextflow processes are
    #   run inside their own directories automatically)
    # Note that we do not need to move the files afterwards with Nextflow, as
    #   the output files will automatically be published in the directory
    #   specified `publishDir`.)
    """
}

process run_multiqc {
    publishDir "${resultsdir}/quality-controls/",
        mode: "copy"

    input:
    path(zips)

    output:
    path("*.html")
    path("multiqc_data/multiqc_general_stats.txt")

    script:
    """
    multiqc -n multiqc.html .

    # Note that we do not need to move/rename any output files.
    # Note that we do not need to remove any of the intermediate directories
    #   that MultiQC creates, as they are published and can easily be cleanup up
    #   by Nextflow or by removing the `work/` directory.
    """
}
