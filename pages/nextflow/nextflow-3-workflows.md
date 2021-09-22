The most important file for any Nextflow workflow is the `main.nf` file, in
which the workflow itself is written. You can view the Nextflow implementation
of the MRSA pipeline in the `tutorials/nextflow/main.nf` file, through which we'll 
go piece by piece in this tutorial. Let's start at the beginning!

These are the first few lines of the file `main.nf`:

```groovy
// Enable DSL2 for a more modular workflow
nextflow.enable.dsl=2

// Define where results should end up
resultsdir = "results/"
```

The first part enables modular designs of Nextflow pipelines (more on this
later), while the second simply defines where the workflow's results should end
up. The next part really defines the workflow itself:

```groovy
workflow {
    // Workflow for generating count data for the MRSA case study

    // Define SRA input data channel
    Channel
        .fromList( params.sra_id_list )
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
```

Let's break it down! The whole paragraph contained in `workflow{}` is the
entirety of the workflow, giving you an overview of it all at a glance. The
first line contained in the triple quotes is just a comment describing the
pipeline, and doesn't actually do anything other than provide a description.

After this comes the definition of the `sra_ids` channel for input data. Line
by line, we first create a new channel using the `Channel` directive, define
that the input should be taken as a list (using `.fromList()`) with the values
stored in the `sra_id_list` parameter, and finally `set` the channel name to
`sra_ids`. The parameters themselves come from the `nextflow.config` file, which
contains this:

```groovy
params {
    sra_id_list = ["SRR935090", "SRR935091", "SRR935092"]
}
```

The next part is the definition of the workflow itself. If you're used to functional
programming you may notice that it looks very much like functions with
arguments, which is the previously mentioned modularity of Nextflow. The first
process is called `get_sra_by_accession` and it takes one input argument; the 
previously defined `sra_ids` input channel. Since we have three separate values 
in this channel, it means that the `get_sra_by_accession` process will be executed 
three times. 

What follows are the remaining processes and their respective input arguments.
You can see that we can take the output from one process and use it as the input
for another. For example, the `run_fastqc` process uses the output from the
`get_sra_by_accession` process, which is done using the `.out` call. Some
processes have no input arguments (such as the `get_genome_fasta` process) while
others use more than one input argument.

The only other thing that is used here is the `collect()` operator, which
collects a channel's content into a single stream. The `run_multiqc` process
uses this operator to collect all the output from the `run_fastqc` process; the
`run_fastqc` process will thus be run three times (one for each sample), while
the `run_multiqc` process will only be run once.

> **Quick recap:** <br>
> In this section we covered:
>
> - Defining channels for input data
> - Defining parameters in a separate configuration file
> - Defining workflows as sets of function-like processes with arguments
