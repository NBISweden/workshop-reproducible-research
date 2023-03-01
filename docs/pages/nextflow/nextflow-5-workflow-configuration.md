We've so far been working with a relatively non-generalised workflow: it's got
hard-coded inputs, paths and genome references. This is perfectly fine for a
project that is purely aimed at getting reproducible results (which is the full
extent of what you want in a lot of cases), but it can be made a lot more
generalisable. Let's go through the MRSA workflow and see what can be improved!

# Configuration basics

One of the things that allow generalisability of Nextflow workflows is
*parameters*, which hold information and values that can be changed directly on
the command-line at the time of execution. One use of parameters in our MRSA
workflow is to remove the hard-coded `results` output directory, for example.
Parameters can be written in the following form:

```nextflow
params {
    parameter_1 = "some/data/path"      // A string parameter
    parameter_2 = 42                    // A value parameter
    parameter_3 = ["a", "b", "c", "d"]  // A list parameter
}
```

You would then refer to these parameters using *e.g.* `params.parameter_1`
anywhere you need to in the workflow. Although parameters can be defined in
`main_mrsa.nf`, it is preferable to define them in a separate *configuration file*. The default name
of this file is `nextflow.config` and if such a file is present it will be used
automatically by Nextflow (to supply a config file with another name use
`nextflow -c <path-to-config-file> run main_mrsa.nf`)

* Create a configuration file and add a parameter for the `results` output
  directory.

* Use your newly created parameter in the `publishDir` directory of a process
  (it'll be in the form of `${params.resultsdir}/some/other/path`, for example).
  Run your workflow to see if it worked.

!!! Tip
    Instead of manually changing all the hard-coded directories in your
    workflow you can use the following little `sed` command, which will do it for
    you: `sed 's/\"results\//\"${params.resultsdir}\//g' main_mrsa.nf > tmp; mv
    tmp main_mrsa.nf`. In case you used a parameter name other than `resultsdir`
    update the command accordingly.

# Command line parameters

Workflow parameters can be assigned on the command-line by executing workflows like
so: `nextflow run main_mrsa.nf --parameter_name 'some_value'`. The workflow parameter `parameter_name`,
is prefixed by a double dash `--` to tell Nextflow this is a parameter to the workflow (a single dash
is a parameter to Nextflow, e.g. `-resume`). The value is also quoted (this is important for parameters
that take file paths as values).

* Run your workflow using the parameter you previously created, but pick
  something other than the default value!

You should now have a new directory containing all the results! This is highly
useful if you want to keep track of separate runs of a workflow with different
software parameters, for example: `nextflow run main.nf --important_param 'value1'
--resultsdir 'value1'`, or simply want to keep the results of separate versions of
the same workflow. You can also change parameters by using the `-params-file` option
or by using another configuration file (and using `-c`), rather than on the command line!

# Configuring inputs

Remember the input for the MRSA workflow, the `ch_sra_ids` channel? This input
is also hard-coded inside the `main_mrsa.nf` file. This could also be made into
a parameter!

* Add another parameter for the input SRA IDs and execute your workflow to check
  that it worked.

Using lists for parameters has its problems though, as you won't be able to
change it on the command line, since the command line doesn't know about Groovy
lists. There are several other ways of specifying inputs in a command
line-friendly way, one of which is to use *sample sheets*. Instead of
specifying samples directly in the command line, you specify a file that lists
the input samples instead; this is the standard that is used in *e.g.* [nf-core
pipelines](https://nf-co.re/). Such a sample sheet for the MRSA workflow might
be stored in *e.g.* `input.csv` and look like this:

```no-highlight
sra_id
SRR935090
SRR935091
SRR935092
```

Reading input from a CSV file can be done by combining the `.splitCsv` channel
factory (splitting the rows to read each entry) and the `.map` operator
(defining which columns should be used). Let's see if we can make it work!

* Create the `input.csv` file with the above shown content.

* Change the definition of the `ch_sra_ids` channel to take the value of a new
  parameter of your choice, defined in the configuration file.

* Add the `.splitCsv(header: true)` operator to the end of the channel
  definition, so that the input is read from the file contents.

* Add the `.map{row -> row.sra_id}` operator, which specifies that each row
  should contain the `sra_id` column, but no other columns.

You should now have a more generalised input to your workflow! Try to run it to
make sure it works - look below if you need some help.

??? example "Click to show the solution"
    ```nextflow
    // Channel definition
    ch_sra_ids = Channel
        .fromPath ( params.sra_ids )
        .splitCsv ( header: true )
        .map      { row -> row.sra_id }

    // Configuration file
    sra_ids = "input.csv"
    ```

By specifying inputs from sample sheets like this we can change inputs
of a workflow execution by creating another sample sheet and specifying
*e.g.*, `--sra_ids input-2.csv` on the command line. This is highly useful when
you want to run a single sample *e.g.*, when testing a workflow, or when you
want to keep track of all the different inputs you've used historically. Sample
sheets are also useful for keeping other metadata, such as custom sample names,
sample groups, location of files, *etc.* For example:

```no-hightlight
sample-1,case,data/sample-1.fastq.gz
sample-2,ctrl,data/sample-2.fastq.gz
sample-3,case,data/sample-3.fastq.gz
sample-4,ctrl,data/sample-4.fastq.gz
```

Here we have not only names and file paths but also to which group each sample
belongs, *i.e.* case or control. Such metadata can be highly useful for more
advanced workflows to use in downstream analyses, such as differential gene
expression! We could create a tuple based on this metadata like so:

```nextflow
ch_input = Channel
    .fromPath("metadata.csv")
    .splitCsv(header: ['id', 'group', 'fastq'])
    .view()
```

> **Input file formatting** <br>
> The input file may also have headers, in which case you can use `header:
> true` to use the column headers defined in the file. You can also read *e.g.*
> tab-separated files by using the `sep` field: `.splitCsv(sep: "\t")`.

# Other configuration scopes

There are lots of things that you might want to add to your configuration, not
just parameters! The workflow *manifest*, for example, which might look like
this:

```nextflow
manifest {
    name        = "My Workflow"
    description = "My workflow, created by me"
    author      = "Me"
    mainScript  = "main.nf"
    version     = "1.0.0"
}
```

* Go ahead and add a workflow manifest to your `nextflow.config` file!

The manifest is useful when you're publishing or sharing the workflow through
*e.g.* GitHub or similar. There are many more such configuration *scopes* that
you might want to use - read more about them [in the documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes).

!!! Success "Quick recap"
    In this section we learnt:

    * How to create parameters in a configuration file
    * How to specify parameters on the command line
    * How to add workflow manifest and other configuration scopes
