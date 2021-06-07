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
    get_genome_fasta()
    index_genome(get_genome_fasta.out)
    align_to_genome(get_sra_by_accession.out, index_genome.out)
    sort_bam(align_to_genome.out)
    get_genome_gff3()
    generate_counts_table(sort_bam.out.collect(), get_genome_gff3.out)
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

process get_genome_fasta {
    publishDir "${resultsdir}/",
        mode: "copy"

    output:
    path("*.fa.gz")

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
    """
}

process index_genome {
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    input:
    path(fasta)

    output:
    path("*.bt2")

    script:
    """
    # Bowtie2 cannot use .gz, so unzip to a temporary file first
    gunzip -c ${fasta} > tempfile
    bowtie2-build tempfile NCTC8325
    """
}

process align_to_genome {
    tag "${sample}"
    publishDir "${resultsdir}/",
        mode: "copy"

    input:
    tuple val(sample), path(fastq)
    path(idx)

    output:
    tuple val(sample), path("*.bam")

    script:
    """
    bowtie2 -x NCTC8325 -U ${fastq} > ${sample}.bam
    """
}

process sort_bam {
    tag "${sample}"
    publishDir "${resultsdir}/",
        mode: "copy"

    input:
    tuple val(sample), path(bam)

    output:
    path("*.sorted.bam")

    script:
    """
    samtools sort ${bam} > ${sample}.sorted.bam
    """
}

process get_genome_gff3 {
    publishDir "${resultsdir}/",
        mode: "copy"

    output:
    path("*.gff3.gz")

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz
    """
}

process generate_counts_table {
    publishDir "${resultsdir}/",
        mode: "copy"

    input:
    path(bam)
    // path(bai)
    path(annotation)

    output:
    path("counts.tsv")

    script:
    """
    # The transcript name is annotated as "Name" in the GFF
    featureCounts -t exon -g Name -a ${annotation} -o counts.tsv ${bam}
    """
}
