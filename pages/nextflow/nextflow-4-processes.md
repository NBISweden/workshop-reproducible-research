Let's look at the first process definition of our `main.nf` Nextflow file:

```groovy
process get_sra_by_accession {
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run
    accession number.
    """
    tag "${sra_id}"
    publishDir "${resultsdir}/data/",
        mode: "copy"

    input:
    val(sra_id)

    output:
    tuple val(sra_id), path("${sra_id}.fastq.gz")

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

The `tag` directive defines what should be shown during the execution of the
process as it is run, *i.e.* the SRA ID in this case. This is useful for viewing
the currently running tag (sample) for any given process at a glance; you can
see exactly how this looks if you run the pipeline one more time with the command
`nextflow run .`.

The `publishDir` directive is used to define where the process' outputs should
be placed. Remember that Nextflow runs each process in its own hidden directory
where only files relevant to its execution exist? That means that we actually
need to move or copy the output files from that hidden directory after a
successful execution of the process to wherever we want them to be stored, which
is what `publishDir` is for: we specify a path (here derived from the
`resultsdir` variable we defined in the beginning of the file) and how the files
should be handled (copied, in this case).

Next comes the `input` directive. We first specify that the input is a value
using `val()`, since it's coming from the `sra_ids` channel we previously
defined, and we name the input variable `sra_id`.

Then comes the `output` directive, which is defined as a `tuple`, *i.e.* having
more than one entry. The output of this process is a combination of the
sample name (from the `sra_id` value variable) and the FASTQ file containing the
reads for that sample. Nextflow will look for files corresponding to the path
defined here and output them to the `publishDir` directory, as well as pass the
entire output tuple to any downstream process that uses them.

The last part of the process is the `script` directive, which works in the same
way as for Snakemake: you either write something in bash itself or call some
external script you've defined elsewhere. Nextflow variables are called using
the syntax `${NF_VARIABLE}`, while bash variables need to be preceded with a
backslash, like so: `\${BASH_VARIABLE}`.

> **Using external scripts** <br>
> While not used in this example, if you have an external script you can put it 
> in a `bin/` directory in the same directory as the `main.nf` file, which will
> automatically give all processes access to that script without you needing to
> specify its absolute path - convenient!

Let's look at the next process, which is the next in line in our workflow
definition:

```groovy
process run_fastqc {
    """
    Run FastQC on a FASTQ file.
    """
    tag "${sample}"
    publishDir "${resultsdir}/qc/",
        mode: "copy",
        saveAs: { filename ->
            filename.indexOf(".zip") > 0 ? \
                "intermediate/${filename}" : "${filename}"
        }

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
two separate outputs? The `saveAs` parameter of the `publishDir` directive
allows us to conditionally publish the output files in different locations, in
this case depending on whether the file ends in `.zip` or not.

> **Syntax explanation** <br>
> The Groovy syntax used here translates to roughly `if filename contains
> ".zip", publish it in the "intermediate/" directory, otherwise publish it in
> the default directory`.

If we move on to the `output` directive, we find another new parameter: `emit` 
is used to name the two different outputs, so that they may be referred to 
separately in downstream processes. For example, the
`run_multiqc` process uses the `run_fastqc.out.zip` as input, *i.e.* only the
compressed files from this process' output.

> **Note** <br>
> Notice that you can both use specific files in output definitions (*e.g*
> `{$sra_id}.fastq.gz`) as well as any file defined by some wildcard, commonly
> using a file extension (*e.g.* `*.zip`). This is useful when you want
> only a specific file or several files.

Remember that this process was called as `run_fastqc(get_sra_by_accession.out)`,
which means that the output tuple from the `get_sra_by_accession` will be the
input for the `run_fastqc` process. This allows us to easily grab both the
sample name and the corresponding FASTQ file without any use of wildcards, which
is how you would do it in Snakemake.

> **Variable names** <br>
> Remember that we can use any name for each input variable within a process,
> regardless of what the variable is called outside the process, just like in
> functional programming. We use the name `sample` here, but the same value was
> called `sra_id` in the previous process.

The script part of this process is slightly changed from its Snakemake
counterpart: firstly, we don't need to use the `-o` flag as all processes in
Nextflow are run in their own directory; and secondly, we don't need to move the
files afterwards, as they are automatically published in the directory specified
by `publishDir`. The same kind of difference holds true for the `run_multiqc`
process as well.

> **Quick recap:** <br>
> In this section we covered:
>
> - The `tag`, `publishDir`, `input`, `output` and `script` process directives
> - How to publish different file types in separate output directories
> - How to give different variable names to separate outputs using `emit`
