Now that we've gone through the specifics of executing workflows in a bit more
detail, let's go through working with processes. While there are numerous
process directives that can be used, we'll go through some of the more commonly
used ones here.

# Tags

Let's look at the command line output we got during the workflow's execution,
which should look something like this:

```bash
N E X T F L O W  ~  version 22.10.6
Launching `./main.nf` [friendly_bhaskara] - revision: b4490b9201
executor >  local (17)
[c9/e5f818] process > DONWLOAD_FASTQ_FILES (SRR935092) [100%] 3 of 3 ✔
[d5/b5f24e] process > RUN_FASTQC (SRR935092)           [100%] 3 of 3 ✔
[91/2cea54] process > RUN_MULTIQC                      [100%] 1 of 1 ✔
[e0/b4fd37] process > GET_GENOME_FASTA                 [100%] 1 of 1 ✔
[87/32ce10] process > INDEX_GENOME                     [100%] 1 of 1 ✔
[56/e9a460] process > ALIGN_TO_GENOME (SRR935092)      [100%] 3 of 3 ✔
[ed/d8c223] process > SORT_BAM (SRR935092)             [100%] 3 of 3 ✔
[e7/4a6bda] process > GET_GENOME_GFF3                  [100%] 1 of 1 ✔
[e9/84f093] process > GENERATE_COUNTS_TABLE            [100%] 1 of 1 ✔
```

Have you noticed that there are SRA IDs after some of the processes? Well, if
you look at which processes show these SRA IDs you might see that it's only
those processes that are executed three times, *i.e.* once per SRA ID. This
doesn't happen automatically, however, and comes from something called *tags*.
Let's look at the `DONWLOAD_FASTQ_FILES` process:

```nextflow
process DONWLOAD_FASTQ_FILES {

    // Download a single-read FASTQ file from the SciLifeLab Figshare remote

    tag "${sra_id}"
    publishDir "results/data",
        mode: "copy"

    input:
    tuple val(sra_id), val(figshare_link)

    output:
    tuple val(sra_id), path("*.fastq.gz")

    script:
    """
    wget ${figshare_link} -O ${sra_id}.fastq.gz
    """
}
```

You can see the `tag` directive at the very top of the process definition. Tags
can be used to *e.g.* show information about the sample currently being analysed
by the process. This is useful both during run-time (allowing you to see which
sample is being processed) but also for debugging or finding problematic samples
in case of errors or odd output. There is, naturally, no need to use tags for
processes which are only run once.

* Comment out (prefix with `//`) the `tag` directive from the
  `DONWLOAD_FASTQ_FILES` process and run the workflow again. What do you see?

Without the `tag` directive you should instead see the numbers 1 through 3,
representing the input files (of which there are three). Nextflow still tells us
that it's working on one of the input files, but it's generally much more useful
to actually see the sample name or ID, rather than just a number.

* Uncomment the `tag` directive before you move on.

# Named outputs

Let's move on to the next process! It looks like this:

```nextflow
process RUN_FASTQC {

    // Run FastQC on a FASTQ file.

    tag "${sample}"
    publishDir "results/",
        mode: "copy"

    input:
    tuple val(sample), path(fastq)

    output:
    path("*.html")
    path("*.zip")

    script:
    """
    fastqc ${fastq} -q
    """
}
```

Here is a process with two output channels! One contains all the `.html` files,
while the other contains all the `.zip` files. How is this handled in the
workflow definition of downstream processes that use the outputs? The
`RUN_MULTIQC` process uses this output, and its part in the workflow definition
looks like this:

```nextflow
RUN_MULTIQC (
    RUN_FASTQC.out[1].collect()
)
```

We already know about `.out` and `.collect()`, but we have something new here:
the `RUN_MULTIQC` process is taking the second channel of the output from the
`RUN_FASTQC` process - `[1]` is the index for the second channel, as Groovy is
zero-based (the first channel is indexed by `[0]`).

This comes with some issues, however. What if we accidentally changed the order
of the outputs in the rule, or added a new one? Using positions like this is
easy to mess up, but there is a better solution: named outputs! This can be
achieved by adding the `emit` option for some or all of the outputs, like so:

```nextflow
output:
path(*.txt), emit: text
```

Instead of referring to the output by its position in an array as before we
refer to the channel with a label we choose (`.out.text`) instead. This benefits
us in that we can infer more information about channel contents called `text`
rather than `[1]`, and it is also allows us to be less error-prone when
rewriting parts of a workflow.

* Your turn! Add named outputs to the `RUN_FASTQC` process and make
  `RUN_MULTIQC` use those outputs. You'll have to change both the output section
  of the `RUN_FASTQC` process, and the workflow definition section for
  `RUN_MULTIQC`. If you need help, see the hint below.

<details>
<summary> Click to show </summary>

```nextflow
// Workflow definition for RUN_MULTIQC
RUN_MULTIQC (
    RUN_FASTQC.out.zip.collect()
)

// Output section of RUN_FASTC
output:
path("*.html"), emit: html
path("*.zip"),  emit: zip
```

</details>

Check if it works by executing the workflow.

# Advanced publishing

So far we've only used the `publishDir` directive in a very simple way:
specifying a directory and the `mode` to use when publishing (to copy the files
rather than symbolically link them). There are more things you can do, however,
especially for processes with more than one output. For example, we can publish
outputs in separate directories, like so:

```nextflow
publishDir "results/tables",
    pattern: "*.tsv",
    mode: "copy"
publishDir "results/logs",
    pattern: "*.log",
    mode: "copy"
```

In this example, `*.tsv` files are copied to the folder `results/tables/`,
while `*.log` files are copied to the folder `results/logs`. The
`publishDir` directive can be used multiple times in a single process, allowing
one to separate output as above, or publish the same output to multiple folders.

* Edit the `RUN_FASTQC` process to place the HTML and compressed files in
  separate directories. Remove the `results` directory and re-run the workflow
  to check that it worked - click below if you're having trouble.

<details>
<summary> Click to show </summary>

```nextflow
process RUN_FASTQC {

    (...)

    publishDir "results/fastqc/html",
        pattern: "*.html",
        mode: "copy"
    publishDir "results/fastqc/zip",
        pattern: "*.zip",
        mode: "copy"

    (...)
}

```

</details>

Note that an output and a *published* output are different things: something can
be an output of a process without being published. In fact, the `RUN_FASTQC`
process is a prime example of this! Think about the compressed output: this
output is only used by the downstream process `RUN_MULTIQC` and is never meant
to be viewed by a human or used by a human in some downstream task not part of
the pipeline itself. We would thus like to keep the compressed files as an
output, but not publish said output. How do we do this? Just remove the
corresponding `publishDir` directive!

The MRSA workflow we've made here was refactored directly from its original
version in the Snakemake tutorial of this course, which means that its output
structure is not fully taking advantage of some of Nextflow's functionality. The
compressed output we've already talked about above is one example.

* See if you can find any other processes in the current implementation of the
  MRSA workflow that you could optimise like this!

Think about whether all processes actually need to have published outputs. Make
sure you test executing the workflow after you've made any changes; click below
if you want a hint.

<details>
<summary> Click to show </summary>

The `GET_GENOME_FASTA` and `GET_GENOME_GFF3` both download reference files which
are only needed by the workflow itself and does not need to be published, the
same goes for the genome index generated by the `INDEX_GENOME` process.

One could argue that neither of the BAM files generated by the `ALIGN_TO_GENOME`
and `SORT_BAM` processes are needed by the user if only the final counts table
is of interest, but BAM files can also be useful for exploring the alignments in
*e.g.* IGV. Both BAMs are, however, definitely not needed: only the sorted one
should be published if one is interested in BAM files.

</details>

# Debugging

It is, sadly, inevitable that we all make mistakes while coding - nobody's
perfect! Nextflow helps you quite a bit when this happens, not just with its
logs but also with informative error messages. Let's introduce an error and look
at what we get:

* Change the final `output` line in the `RUN_MULTIQC` process to the following
  and re-run the workflow: `path("multiqc_general_stats.csv")` - notice the
  usage of `.csv` rather than `.txt` as before.

We got an error! We get a number of things, actually, including (in order from
the top) the name of the process that gave the error, the likely cause, the
command that was executed, along with its exit status, output, error and the
work directory that the task was run in. Let's focus on the `Caused by:` part at
the top, which should look something like this:

```no-highlight
Caused by:
  Missing output file(s) `multiqc_general_stats.csv` expected by process `RUN_MULTIQC`
```

We can also see that the command's exit status is `0`, which means that the
command was successful; any exit status other than `0` means there was an error
of some kind. We can thus infer that the command (1) worked, (2) failed to give
us the output expected by Nextflow. Thankfully, Nextflow graciously prints the
work directory for us so that we may check out what happened in more detail.

* Copy the working directory path, `cd` into it and list its contents using
  `ls`.

You might already have spotted the error in the message above; the error we
introduced here was that the expected output file has a `.csv` extension,
rather than the correct `.txt`. Nextflow is expecting the `.csv` output, but
the process `script` directive is (correctly) giving us the `.txt` file, which
we can see inside the process' work directory.

* Go back to the root directory, revert the error you introduced and re-run the
  workflow to make sure it works again.

This might have seemed like a trivial error, but a lot of errors in Nextflow can
be solved in the same manner, *i.e.* by just following the debugging output
reported by Nextflow and inspecting the specific subdirectory in question.

> **A note about Bash** <br>
> If you are using Bash variables inside the `script` directive you have to be
> careful to prepend it with a backslash, *e.g.* `\${BASH_VARIABLE}`. This is
> because the dollar-sign is used by Nextflow, so you have to tell Nextflow
> explicitly when you're using a Bash variable. This is a common source of
> errors when using Bash variables, so keeping it in mind can save you some
> debugging time!

> **Quick recap** <br>
> In this section we've learnt:
>
> * How to use the `tag` directive
> * How to use named output with `emit`
> * How to publish outputs into different directories
> * How to debug errors and mistakes
