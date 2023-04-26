It's time to start working with a more realistic workflow using the MRSA case
study of this course! We've created a bare-bones version of this pipeline for
you, but we'll work our way through it as we go along and learn more about
Nextflow's features and functionality. The MRSA workflow looks like this:

```nextflow
workflow {

    // Workflow for generating count data for the MRSA case study

    // Get input files from a samplesheet
    ch_input = Channel
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
}
```

The workflow has one input channel named `ch_input`, which reads input from the
`samplesheet.csv` file. We then define the processes to be executed by this
workflow, nine in total. The first process (`DOWNLOAD_FASTQ_FILES`) takes the
`ch_input` channel as input, while the rest of the processes takes the output of
previous processes as input. Before we go into more detail regarding the
ins-and-outs of this workflow, let's start with some specifics of how workflows
are executed and what you can get from them.

# Reports and visualisations

Let's start with running the workflow plus getting some reports and
visualisation while we're at it!

* Run the workflow using the following command: `nextflow run main_mrsa.nf
  -with-report report.html -with-timeline timeline.html -with-dag dag.png`.

After successful executing, you will find three more files in your current
directory: `report.html`, `timeline.html` and `dag.png`. The first file contains
a workflow report, which includes various information regarding execution such
as runtime, resource usage and details about the different processes. The second
file contains a timeline for how long each individual process took to execute,
while the last contains a visualisation of the workflow itself.

Take a few minutes to browse these files for yourself. When running a workflow
you can of course choose which of these additional files you want to include by
picking which ones are important or interesting to you - or don’t include any!

# Logs

Nextflow keeps a log of all the workflows that have been executed. Let's check
it out!

* Type `nextflow log` to get a list of all the executions.

Here we get information about when the workflow was executed, how long it ran,
its run name, whether it succeeded or not and what command was used to run it.
You can also use `nextflow log <run name>` to show each task's directory that
was executed for that run. You can also supply the `-f` (or `-fields`) flag
along with additional fields to show.

* Run `nextflow log <run name> -f hash,name,exit,status`

This shows us not only the beginning of each task's working directory, but also
its name, exit code and status (*i.e.* if it completed successfully or failed in
some manner).

> **Listing fields** <br>
> If you want to see a complete list of all the fields you might explore using
> the log, just type `nextflow log -l` or `nextflow log -list-fields`. This is
> highly useful for debugging when there's some specific information about a run
> you're particularly interested in!

We can also get even more detailed information about the latest
run by looking into the `.nextflow.log` file!

* Look into the latest log by typing `less .nextflow.log`.

You'll be greeted by a wealth of debugging information, which may even seem a
bit overkill at this point! This level of detail is, however, quite useful both
as a history of what you've attempted and as an additional help when you run
into errors! Also, it helps with advanced debugging - which we'll get into
later.

# Re-running workflows

Something you often want to do in Nextflow (or any WfMS for that matter) is to
re-run the workflow when you changed some input files or some of the code for
its analyses, but you don't want to re-run the entire workflow from start to
finish. Let’s find out how this works in Nextflow!

* Run the same `nextflow run main_mrsa.nf` command again.

You'll notice that Nextflow actually re-ran the entire workflow from scratch,
even though we didn't change anything. This is the default behaviour of
Nextflow.

* Let’s try that again: `nextflow run main_mrsa.nf -resume` instead.

Now you can see that Nextflow didn't actually re-run anything. The `-resume`
flag instructed Nextflow to use the cached results from the previous run!

Nextflow automatically keeps track of not only changes to input files, but also
changes to code, process definitions and scripts. You can thus change anything
relating to your workflow and just re-run with the `-resume` flag and be sure
that only processes relevant to your changes are executed again!

* Use `tree work/` to list the contents of the work directory.

Because Nextflow keeps track of all the runs, we've now got two sets of files
in the work directory. One set from the first run, and another from the second
run. This can take up valuable space, so let's clean that up.

* Use `nextflow clean -n -before <run_name>` to show which work directories will
  be cleaned up (use `nextflow log` to find the run name if you don't remember
  it). Then delete those directories by changing `-n` (dry-run) to `-f` (force).

Nextflow's `clean` subcommand can be used to clean up failed tasks and unused
processes. Here we used the `-before` flag, meaning that any runs before the
specified run are removed; use `nextflow help clean` to see other options for
cleaning. This is the preferred way to clean up the working directory.

* Remove the `results` directory and re-run the workflow again using the
  `-resume` flag.

We removed all the results we used before, but we still managed to resume the
workflow and use its cache - how come? Remember that Nextflow uses the `work`
directory to run all of its tasks, while the `results` directory is just where
we have chosen to publish our outputs. We can thus delete the `results`
directory as often as we like (a necessity when output filenames are changed)
and still get everything back without having to re-run anything. If we were to
delete the `work` directory, however...

* Delete the `work` directory and re-run the workflow using the `-resume` flag.

There is no longer any cache for Nextflow to use, so it re-runs from the start!
This is good to keep in mind: you can always delete the output directories of
your workflow, but if you mess with `work` you'll lose, well... work!

> **Quick recap** <br>
> In this section we've learnt:
>
> * How to get automatic reports and visualisations
> * How to check the Nextflow logs
> * How to re-run workflows
> * How to clean the Nextflow cache
