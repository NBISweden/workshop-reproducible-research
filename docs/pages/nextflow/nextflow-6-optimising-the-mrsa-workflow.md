We just added several parameters and configurations to our MRSA workflow, but
we didn't do anything about the reference genomes: those are still hard-coded.
How come? Well, the current MRSA workflow is, in fact, not very well-optimised
for Nextflow at all, being a refactor from the Snakemake tutorial of this
course.

All of the processes are basically unchanged, excluding some minor alterations.
For example, the `run_fastqc` rule in Snakemake used the `-o` flag to specify
that the results should be in the current directory, followed by moving the
output files to their respective output directory. The first part is not needed
in Nextflow (as everything is run in its own subdirectory), and the second part
is done by the `publishDir` directive. These are just minor alterations,
though, but we can do much more if we fully utilise Nextflow's features!

# Remote files

One of these features is the ability to automatically download remote files,
without needing to explicitly do so! The `path` input type can handle either
file paths (like we've done so far) or a URI-supported protocol (such as
`http://`, `s3://`, `ftp://`, *etc.*). This would be highly useful for *e.g.*
the `GET_GENOME_FASTA` process - in fact, we don't even need that process at
all! All we need to do is to change the input to the `INDEX_GENOME` and
`ALIGN_TO_GENOME` processes.

* Create a new input channel using the `fromPath()` channel factory and the
  absolute path to the genome FASTA.

* Make the `INDEX_GENOME` process use that input channel instead of the
  previously used output of the `GET_GENOME_FASTA` process.

* Remove the `GET_GENOME_FASTA` process, as it is not needed anymore.

Re-run the workflow to see if it worked. Check the code below for an example if
you're stuck:

??? example "Click to show the solution"
    ```nextflow
    # Channel creation

    ch_genome_fasta = Channel.fromPath( "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz" )

    # Workflow definition
    INDEX_GENOME (
        ch_genome_fasta
    )
    ```

If we want to get detailed we can also change the hard-coded "NCT8325"
naming in *e.g.* the `INDEX_GENOME` process and put that in another parameter,
or grab the `baseName()` from the channel and make a `[prefix, file]` tuple
using the `map{}` operator like we did previously; check below if you're
curious of how this could be done.

??? example "Click to show the solution"
    ```nextflow
    // Channel definition
    ch_genome_fasta = Channel
        .fromPath( "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna/Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz" )
        .map{ file -> tuple(file.getBaseName(), file) }

    // INDEX_GENOME process definition
    process INDEX_GENOME {

        publishDir "results/intermediate/",
            mode: "copy"

        input:
        tuple val(fasta_name), path(fasta)

        output:
        path("*.b2t"), emit: index

        script:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c ${fasta} > tempfile
        bowtie2-build tempfile ${fasta_name}
        """
    }
    ```

# Subworkflows

The DSL2 allows highly modular workflow design, where a workflow may contain
multiple *subworkflows*. A subworkflow is just like a normal workflow, but it
can be called inside other workflows, similar to a process. There is thus no
special difference between a subworkflow and a workflow; the only difference is
how you use them in practice. Let's take a look at a toy example:

```nextflow
workflow {
    ch_input = Channel.fromPath( params.input )
    SUBWORKFLOW( ch_input )
}

workflow SUBWORKFLOW {

    take:
    input_file

    main:
    ALIGN_READS( input_file )

    emit:
    bam = ALIGN_READS.out.bam
}
```

Here we have an unnamed, main workflow like before, plus a named subworkflow.
A workflow can have inputs specified by the `take` directive, which is
the equivalent of process `input` for workflows. The `main` part is
the workflow body, which contains how to run which processes in which order.
The last part, `emit`, also works the same as for processes, in that we name
the different outputs of the workflow so that we may use them in other
workflows or processes. Nextflow will run the unnamed workflow by default,
unless the `-entry` flag is specified, like so:

```bash
nextflow run main.nf -entry SUBWORKFLOW
```

This will run the workflow named `SUBWORKFLOW`, but nothing else. You can also
store subworkflows in separate files, so that everything doesn't have to be
crammed into a single `main.nf` file. A subworkflow named `SUBWORKFLOW`
contained in the file `subworkflow.nf` can be loaded into a `main.nf` file like
so:

```nextflow
include { SUBWORKFLOW } from "./subworkflow.nf"
```

If you have a complex workflow with several subworkflows you might thus store
them in a separate directory, *e.g.* `subworkflows`. This allows you to have
fine-grained control over the general architecture of your Nextflow workflows,
organising them in a manner that is easy to code and maintain. A `process`
can also be treated in the same manner, and defined separately in another file.

* Now it's your turn! Separate the `RUN_FASTQC` and `RUN_MULTIQC` processes out
  of the main workflow and into a subworkflow. Check below if you're having
  trouble.

??? example "Click to show the solution"
    ```nextflow
    // In the main workflow:
    RUN_QC (
        GET_SRA_BY_ACCESSION.out
    )

    // A new subworkflow
    workflow RUN_QC {

        take:
        fastq

        main:
        RUN_FASTQC (
            fastq
        )
        RUN_MULTIQC (
            RUN_FASTQC.out.zip.collect()
        )
    }
    ```

If you want to challenge yourself, try to do the same with the `INDEX_GENOME`,
`ALIGN_TO_GENOME` and `SORT_BAM` processes! Be careful of where you get your
inputs and outputs.

!!! Success "Quick recap"
    In this section we learnt:

    * How to automatically download remote files
    * How to create and work with subworkflows
