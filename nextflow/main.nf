#!/usr/bin/env nextflow

// Enable DSL2 for a more modular workflow
nextflow.enable.dsl=2

// Define where results should end up
resultsdir = "results/"

workflow {
    """ Workflow for generating count data for the MRSA case study """
    // Here we define something that Nextflow denotes a "channel", which is a
    // stream of data. These are used to define the input to be passed through
    // the workflow - meaning SRA IDs defined as coming from the list in the
    // config file (using `fromList`), in this particular case. The `set{...}`
    // part sets the name for the channel, which is then used as input in the
    // first process (see below).
    Channel
        .fromList( params.sra_id_list )
        .set{ sra_ids }

    // Here we define the workflow itself, which in Nextflow is a sequence of
    // function-like calls: each process (which is the equivalent of Snakemake's
    // rules) is called like a function, i.e. with input arguments. Each
    // process' output can be accessed by appending `.out` to it, which is how
    // downstream processes can use the output of those upstream.
    get_sra_by_accession(sra_ids)
    run_fastqc(get_sra_by_accession.out)
    // The `collect()` function is almost equivalent to Snakemake's `expand()`,
    // in that it collects all the outputs from a particular channel into one:
    // in this way, you can make a process run on ALL of the outputs from some
    // other process in a single execution, rather than once sample at a time.
    run_multiqc(run_fastqc.out.zip.collect())
    get_genome_fasta()
    index_genome(get_genome_fasta.out)
    align_to_genome(get_sra_by_accession.out, index_genome.out)
    sort_bam(align_to_genome.out)
    get_genome_gff3()
    generate_counts_table(sort_bam.out.collect(), get_genome_gff3.out)
}

process get_sra_by_accession {
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run
    accession number.
    """
    // The `tag` parameter is what will be shown for each execution of this
    // particular process, e.g. once for each SRA ID in this case.
    tag "${sra_id}"
    // The `publishDir` parameter is where the output files will be "published"
    // once the process has run, i.e. where they will be copied to once the
    // process has been successful. We could omit the publishing of this
    // process' results as it is only data downloaded from the internet, and is
    // not directly useful for the user to have access to.
    publishDir "${resultsdir}/data/",
        mode: "copy"

    // Notice that the input is an object rather than an explicit string, since
    // we have given the inputs above in the `workflow {...}` specification.
    input:
    val(sra_id)

    // Notice how we're actually giving two thing as one output channel: both
    // the SRA ID (a value) and the downloaded FASTQ file (a file path). This
    // makes it easy to keep a file and its associated sample name together,
    // without having to do any string/filename manipulation.
    output:
    tuple val(sra_id), path("${sra_id}.fastq.gz")

    script:
    // Rather than saying e.g. {input.sra_id} as in Snakemake we only specify
    // the variable name in Nextflow.
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
    """
    Run FastQC on a FASTQ file.
    """
    tag "${sample}"
    publishDir "${resultsdir}/qc/",
        mode: "copy",
        // Here we specify that anything ending in `.zip` should be put into the
        // `intermediate/` directory rathern in the default publishDir. One
        // could argue, however, that since the compressed output of FastQC are
        // just intermediate files for MultiQC we do not need to publish them.
        saveAs: { filename ->
            filename.indexOf(".zip") > 0 ? \
                "intermediate/${filename}" : "${filename}"
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
    // name the output data streams and "emit" them individually, which can then
    // be used by other processes specifically.
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:
    // Note that We don't need to specify `-o .` as when using Snakemake, as all
    // Nextflow processes are run inside their own directories automatically).
    // Also note that we do not need to move the files afterwards with Nextflow,
    // as the output files will automatically be published in the directory
    // specified `publishDir`.)
    """
    # Run FastQC
    fastqc ${fastq} -q
    """
}

process run_multiqc {
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
    publishDir "${resultsdir}/qc/",
        mode: "copy"

    input:
    path(zips)

    output:
    path("*.html")
    path("multiqc_data/multiqc_general_stats.txt")

    script:
    // Note that we do not need to move/rename any output files.
    // Also note that we do not need to remove any of the intermediate
    // directories that MultiQC creates, as they are not published and can
    // easily be cleaned up by Nextflow or by removing the `work/` directory.
    """
    multiqc -n multiqc.html .
    """
}

process get_genome_fasta {
    """
    Retrieve the sequence in fasta format for a genome.
    """
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    output:
    path("*.fa.gz")

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
    """
}

process get_genome_gff3 {
    """
    Retrieve annotation in gff3 format for a genome.
    """
    publishDir "${resultsdir}/idx/",
        mode: "copy"

    output:
    path("*.gff3.gz")

    script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz
    """
}

process index_genome {
    """
    Index a genome using Bowtie 2.
    """
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
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    tag "${sample}"

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
    """
    Sort a bam file.
    """
    tag "${sample}"
    publishDir "${resultsdir}/bam/",
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

process generate_counts_table {
    """
    Generate a count table using featureCounts.
    """
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
