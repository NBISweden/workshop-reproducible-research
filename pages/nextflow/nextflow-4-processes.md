Let's look at the first process definition of our `main.nf` Nextflow file:

```groovy
process GET_SRA_BY_ACCESSION {

    // Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run
    // accession number.

    tag "${sra_id}"
    publishDir "${resultsdir}/data/raw_internal",
        mode: "copy"

    input:
    val(sra_id)

    output:
    tuple val(sra_id), path("${sra_id}.fastq.gz"), emit: sra_data

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
```

The overall structure of a Nextflow process is actually quite similar to
that of a Snakemake rule: we have `input`, `output` and `script` sections,
with the additions of `tag` and `publishDir` directives. Let's go through
them from top to bottom.

### Tags

The `tag` directive defines what should be shown during the execution of the
process as it is run, *i.e.* the SRA ID in this case. This is useful for
viewing the currently running tag (*i.e.* sample) for any given process at a
glance.

* Re-run the workflow without using the cache with `nextflow run main.nf` and
  pay extra attention to how the tags for the processes change.

### Publishing results

The `publishDir` directive is used to define where the process' outputs should
be placed, which is called "publishing" in this context. Remember that Nextflow
runs each process in its own directory where only files relevant to its
execution exist? That means that we actually need to move or copy the output
files from that directory after a successful execution of the process to
wherever we want them to be stored, which is what `publishDir` is for: we
specify a path (here derived from the `resultsdir` variable we defined in the
beginning of the file) and how the files should be handled (copied, in this
case).

* Edit the `GET_SRA_BY_ACCESSION` process so that its output is published in the
  `${resultsdir}/data` instead and re-run the workflow using the `-resume` flag.

Notice that Nextflow didn't actually re-run a single process, even though we
changed the `publishDir` directive? This is because we didn't actually change
any functional code, just where the results should be published - Nextflow keeps
track of this difference, so if we want to change to structure of our results
directory that's easily done without having to re-run anything! This is the same
as for when we changed the `resultsdir` variable in the previous tutorial.

### Inputs

Next comes the `input` directive. We first specify that the input is a value
using `val()`, since it's coming from the `ch_sra_ids` channel we previously
defined, and we name the input variable `sra_id`. The input for this particular
process is not that complex, but we'll look at something more complicated in a
little while.

### Outputs

Then comes the `output` directive, which is defined as a `tuple`, *i.e.* having
more than one entry. The output of this process is a combination of the
sample name (from the `sra_id` value variable) and the FASTQ file containing the
reads for that sample. Nextflow will look for files corresponding to the path
defined here and output them to the `publishDir` directory, as well as pass the
entire output tuple to any downstream process that uses them.

The `emit` directive is used to name the output for use in downstream processes:
for example, another process might be called in the workflow definition like
this: `DOWNSTREAM_PROCESS(GET_SRA_BY_ACCESSION.out.sra_data)`, which will then
use the specific `sra_data` output of the `GET_SRA_BY_ACCESSION` process as
input.

> **Note** <br>
> If a process only has a single output you can skip using the `emit` directive,
> and just use `.out` for any downstream process using that output instead of
> `.out.<emitted output>`. We have, however, opted to use the `emit` directive
> for all processes and outputs in this material for clarity and consistency.

### Code and scripts

The last part of the process is the `script` directive, which works in the same
way as for Snakemake: you either write something in bash itself or call some
external script you've defined elsewhere. When you contain some text in triple
quotes inside the `script` directive you will run whatever code you write using
bash, just like in Snakemake.

Using Nextflow variables inside its Bash context works the same way as for
Snakemake: `$nf_variable` (or with squiggly brackets: `${nf_variable}`). This
clashes with Bash variables, however, which works exactly the same way (*e.g.*
`$bash_variable` or `${bash_variable}`). This means that when we refer to a Bash
variable in the Nextflow `script` bash context, we need to prepend the dollar
sign with a backslash, like so: `\${bash_variable}`. In the
`GET_SRA_BY_ACCESSION` process above we are using `${sra_id}`, which tells us
that it is a Nextflow variable, which comes from the process input in this case.

> **Using external scripts** <br>
> While not used in this example, if you have an external script you can put it
> in a `bin/` directory in the same directory as the `main.nf` file, which will
> automatically give all processes access to that script without you needing to
> specify its absolute path - convenient!

### A more complicated process

Let's look at the next process, which is the next in line in our workflow
definition:

```groovy
process RUN_FASTQC {

    // Run FastQC on a FASTQ file.

    tag "${sample}"
    publishDir "${resultsdir}/results/",
        pattern: "*.html",
        mode: "copy"
    publishDir "${resultsdir}/intermediate",
        pattern: "*.zip",
        mode: "copy"

    input:
    tuple val(sample), path(fastq)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:
    """
    # Run FastQC
    fastqc ${fastq} -q
    """
}
```

Let's start with the extended `publishDir` directive. Notice how the process has
two separate outputs, as well as two separate `publishDir` directives? We use
`pattern` as part of the `publishDir` directives to put the HTML and ZIP files
in separate output directives. This is useful when you want to separate the
final location of your different outputs. Another useful instance for this is
when you have intermediate data that you don't necessarily want to publish: one
could argue that this is the case here, where the `.zip`-files are not really
something the end-user wants to look at, but they are still required as output
for the process.

* Remove the second `publishDir` directive, delete the results directory (`rm -r
  results`) and re-run the workflow using `nextflow run main.nf -resume`.

You should no longer see the compressed ZIP-files in the `results/intermediate`
directory, but they are still available for downstream processes that use them.

If we move on to the `output` directive, we find that we have two separate
outputs, each with their own `emit` directive, so that they may be referred to
separately in downstream processes. For example, the `RUN_MULTIQC` process uses
the `RUN_FASTQC.out.zip` as input, *i.e.* only the compressed files from this
process' output.

> **Note** <br>
> Notice that you can both use specific files in output definitions (*e.g*
> `${sra_id}.fastq.gz`) as well as any file defined by some wildcard, commonly
> using a file extension (*e.g.* `*.zip`). This is useful when you want
> only a specific file or several files.

Remember that this process was called as `RUN_FASTQC(GET_SRA_BY_ACCESSION.out.sra_data)`,
which means that the output tuple from the `GET_SRA_BY_ACCESSION` will be the
input for the `RUN_FASTQC` process. This allows us to grab both the sample name
and the corresponding FASTQ file without any use of wildcards, which is how you
would do it in Snakemake.

> **Variable names** <br>
> Remember that we can use any name for each input variable within a process,
> regardless of what the variable is called outside the process, just like in
> functional programming. We use the name `sample` here, but the same value was
> called `sra_id` in the previous process.

The script part of this process is slightly changed from its Snakemake
counterpart: firstly, we don't need to use the `-o` flag as all processes in
Nextflow are run in their own directory; and secondly, we don't need to move the
files afterwards, as they are automatically published in the directory specified
by `publishDir`.

> **Quick recap:** <br>
> In this section we covered:
>
> - The `tag`, `publishDir`, `input`, `output` and `script` process directives
> - How to give different variable names to separate outputs using `emit`
> - How to publish different output files in separate directories
