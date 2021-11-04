The most important file for any Nextflow workflow is the `main.nf` file, in
which the workflow itself is written. You can view the Nextflow implementation
of the MRSA pipeline in the `tutorials/nextflow/main.nf` file, through which
we'll go piece by piece in this tutorial. Let's start at the beginning!

These are the first few lines of the file `main.nf`:

```groovy
// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2

// Define where results should end up
resultsdir = "results/"
```

The first part enables more powerful functionality of Nextflow, including
modular designs of pipelines (more on this later), while the second simply
defines where the workflow's results should end up. Double slashes (`//`) are
used as comments in Nextflow. Let's see how we can change the final output
directory of Nextflow!

* Change the second line above to read `resultsdir = "new-results/"` instead and
  re-run the workflow with `nextflow run main.nf -resume`.

You should now see that we have a `new-results` directory containing all the
workflow's results, identical to the previous ones, and Nextflow didn't even
have to re-run anything! That's because Nextflow runs each process in its own
directory inside `work/`, as previously mentioned, so all it needs to do is to
put the final results in the `new-results/` directory instead. Before reading
on, change the `resultsdir` parameter back to `results/` and remove the
`new-results/` directory, so that we're back to how things were just a moment
ago.

The next part defines the workflow itself:

```groovy
workflow {

    // Workflow for generating count data for the MRSA case study

    // Define SRA input data channel
    Channel
        .fromList ( params.sra_id_list )
        .set      { ch_sra_ids }

    // Define the workflow
    GET_SRA_BY_ACCESSION (
        ch_sra_ids
    )
    RUN_FASTQC (
        GET_SRA_BY_ACCESSION.out.sra_data
    )
    RUN_MULTIQC (
        RUN_FASTQC.out.zip.collect()
    )
    GET_GENOME_FASTA ()
    INDEX_GENOME (
        GET_GENOME_FASTA.out.fasta
    )
    ALIGN_TO_GENOME (
        GET_SRA_BY_ACCESSION.out.sra_data,
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
```

Let's break it down! The whole paragraph contained in `workflow{}` is the
entirety of the workflow, giving you an overview of it all at a glance.

The first part is the definition of the `ch_sra_ids` channel for input data.
Line by line, we first create a new channel using the `Channel` directive,
define that the input should be taken as a list (using `.fromList()`) with the
values stored in the `sra_id_list` parameter, and finally `set{}` the channel
name to `ch_sra_ids`. The parameters themselves come from the `nextflow.config`
file, which contains this:

```groovy
params {
    sra_id_list = ["SRR935090", "SRR935091", "SRR935092"]
}
```

> **Naming channels** <br>
> Notice that we prepend our channel name with `ch_`, which is only done for
> readability and ease of development; there is nothing at all that forces you
> to follow this convention if you don't want to, but it's nice to be able to
> see at a glance which variables are channels and which ones are not.

The next part is the definition of the workflow itself. If you're used to
functional programming you may notice that it looks very much like functions
with arguments, which is the previously mentioned modularity of Nextflow. The
first process is called `GET_SRA_BY_ACCESSION` and it takes one input argument:
the previously defined `ch_sra_ids` input channel. Since we have three separate
values in this channel, it means that the `GET_SRA_BY_ACCESSION` process will
be executed three times.

> **Nextflow and whitespace** <br>
> You might have noticed that whitespace is used liberally in a lot of the
> workflow definition above. Nextflow doesn't care at all about whitespace, so
> go ahead and use it in whatever manner you think is easiest to read and work
> with! The style chosen here (including the choice to use UPPERCASE for all
> processes) is based on [nf-core](https://nf-co.re://nf-co.re/).

What follows are the remaining processes and their respective input arguments.
You can see that we can take the output from one process and use it as the
input for another. For example, the `RUN_FASTQC` process uses the named
`sra_data` output from the `GET_SRA_BY_ACCESSION` process, which is done using
the `.out` call. Some processes have no input arguments (such as the
`GET_GENOME_FASTA` process) while others use more than one input argument,
which is defined in each process (just like how functions work).

The only other thing that is used here is the `collect()` operator, which
collects a channel's content into a single stream. The `RUN_MULTIQC` process
uses this operator to collect all the output from the `RUN_FASTQC` process; the
`RUN_FASTQC` process will thus be run three times (one for each sample), while
the `RUN_MULTIQC` process will only be run once.

* Let's write a workflow of our own, using the already existing processes! Copy
  the following paragraph into the bottom of the `main.nf` file:

```groovy
workflow MULTIQC {

    // Workflow for running FastQC and MultiQC on MRSA data

    // Define SRA input data channel
    Channel
        .fromList ( params.sra_id_list )
        .set      { ch_sra_ids }

    // Define the workflow up until the MultiQC process
    GET_SRA_BY_ACCESSION (
        ch_sra_ids
    )
    RUN_FASTQC (
        GET_SRA_BY_ACCESSION.out.sra_data
    )
    RUN_MULTIQC (
        RUN_FASTQC.out.zip.collect()
    )
}
```

Notice that it's exactly the same as the whole workflow above, except we only
kept the processes up until running MultiQC. We also named the workflow
`MULTIQC`, instead of having no name as above. The unnamed workflow is
considered the "main" workflow and is executed by default, but we can run any
other workflow using the `-entry <WORKFLOW>` flag.

* Delete the previous results with `rm -r results/` and then run the new
  workflow by executing `nextflow run main.nf -entry MULTIQC`.

You should now see that only the processes defined in our MULTIQC workflow were
executed. You may thus use this kind of structure to define sub-workflows that
run only part of your analyses and subsequently string them together in whatever
order you desire.

> **Quick recap:** <br>
> In this section we covered:
>
> - Defining channels for input data
> - Defining parameters in a separate configuration file
> - Defining workflows as groups of function-like processes with arguments
> - Creating and executing sub-workflows
