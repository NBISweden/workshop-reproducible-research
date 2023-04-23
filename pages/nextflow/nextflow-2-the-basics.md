We'll start by creating a very simple workflow from scratch, to show how
Nextflow works: it will take two input files and convert them to UPPERCASE
letters.

* Start by running the following commands:

```bash
touch main.nf
echo "This is a.txt" > a.txt
echo "This is b.txt" > b.txt
```

Open the `main.nf` file with an editor of your choice. This is the main workflow
file used in Nextflow, where workflows and their processes are defined.

* Copy the following code into your `main.nf` file:

```nextflow
nextflow.enable.dsl = 2

// Workflow definition
workflow {
    // Define input files
    ch_input = Channel.fromPath( "a.txt" )

    // Run workflow
    CONVERT_TO_UPPER_CASE( ch_input )
}

// Process definition
process CONVERT_TO_UPPER_CASE {
    publishDir "results/",
        mode: "copy"

    input:
    path(file)

    output:
    path("a.upper.txt")

    script:
    """
    tr [a-z] [A-Z] < ${file} > a.upper.txt
    """
}
```

> **DSL2** <br>
> The `nextflow.enable.dsl = 2` line is only required for earlier versions of
> Nextflow, so depending on which version you have installed it might be
> redundant; all it does is enable the latest functionality of modular Nextflow
> design.

Here we have two separate parts. The first is the *workflow definition*, while
the last is a *process*. Let's go through them both in more detail!

> **Nextflow comments** <br>
> Double-slashes (`//`) are used for comments in Nextflow.

> **Nextflow and whitespace** <br>
> Nextflow is not indentation-sensitive. In fact, Nextflow doesn't care at all
> about whitespace, so go ahead and use it in whatever manner you think is
> easiest to read and work with! Do keep in mind that indentations and other
> types of whitespace *does* improve readability, so it's generally not a good
> idea to forego it entirely, even though you can.

# Workflow definitions

```nextflow
workflow {
    // Define input files
    ch_input = Channel.fromPath( "a.txt" )

    // Run workflow
    CONVERT_TO_UPPER_CASE( ch_input )
}
```

The workflow definition here has two parts, each doing an important job for any
Nextflow workflow. The first part defines a *channel*, which is an asynchronous
first-in-first-out stream of data that connect a workflow's various inputs and
outputs. In simpler terms, channels contain the data that you want to process
with the workflow and can be passed between the various parts of the workflow.

Channels can be created in various different ways using *channel factories*,
depending on what type data you want to put into them and where this data is
stored. In this particular case we define our `ch_input` channel using the
`.fromPath` channel factory, which takes a file path as input - here we use the
`a.txt` file. You can thus read `ch_input = Channel.fromPath("a.txt")` as
*"create the channel `ch_input` and send the file `a.txt` into it"*.

> **Naming channels** <br>
> A channel can be named anything you like, but it is good practice to prepend
> them with `ch_`, as that makes it clear which variables are channels and which
> are just normal variables.

How do we use these channels then? Channels pass data to and from processes
through our workflow. By providing channels as arguments to processes, we
describe how we want data to flow. This is exactly what we do in the second
part: we call our `CONVERT_TO_UPPER_CASE` process with the `ch_input` as input
argument - this is very similar to functional programming.

This is our entire workflow, for now: the creation of a channel followed by
using the contents of that channel as input to a single process. Let's look at
how processes themselves are defined!

# Process definitions

```nextflow
process CONVERT_TO_UPPER_CASE {
    publishDir "results/",
        mode: "copy"

    input:
    path(file)

    output:
    path("a.upper.txt")

    script:
    """
    tr [a-z] [A-Z] < ${file} > a.upper.txt
    """
}
```

Looking at the process in the code above, we can see several parts. The process
block starts with its name, in this case `CONVERT_TO_UPPER_CASE`, followed by
several sections, or *directives* as Nextflow calls them: `publishDir`, `input`,
`output` and `script`.

> **Naming processes** <br>
> A process can be named using any case, but a commonly used convention is to use 
> UPPERCASE letters for processes to visually distinguish them in the workflow. 
> You do not have to follow this if you don't want to, but we do so here.

Let's start with the first directive: `publishDir`. This tells Nextflow where
the output of the process should be placed when it is finished. Setting `mode`
to `"copy"` just means that we want to copy the output files to the publishing
directory, rather than using a symbolic link (which is the default).

The `input` and `output` directives describe the data expected to come through
this specific process. Each line of `input` describes the data expected for each
process argument, in the order used in the workflow. In this case,
`CONVERT_TO_UPPER_CASE` expects a single channel (one line of input), and
expects the data to be filenames ( *i.e.* of type `path`).

Notice that there is a difference between how the inputs and outputs are
declared? The `output` is an explicit string (*i.e.* surrounded by quotes),
while the input is a variable named `file`. This means inputs can be referenced
in the process without naming the data explicitly, unlike the output where the
name needs to be explicit. We'll get back to exactly how this works in just a
moment.

# Executing workflows

Let's try running the workflow we just created!

* Type the following in your terminal:

```bash
nextflow run main.nf
```

This will make Nextflow run the workflow specified in your `main.nf` file. You
should see something along these lines:

```no-highlight
N E X T F L O W  ~  version 22.10.6
Launching `./main.nf` [mad_legentil] - revision: 87f0c253ed
executor >  local (1)
[32/9124a1] process > CONVERT_TO_UPPER_CASE (1) [100%] 1 of 1 ✔
```

The first few lines are information about this particular run, including the
Nextflow version used, which workflow definition file was used, a randomly
generated run name (an adjective and a scientist), the revision ID as well
as where the processes were executed (locally, in this case, as opposed to
*e.g.* SLURM or AWS).

What follows next is a list of all the various processes for this particular
workflow. The order does not necessarily reflect the order of execution
(depending on each process’ input and output dependencies), but they are in the
order they were defined in the workflow file - there's only the one process
here, of course. The first part (*e.g.* `[32/9124a1]`) is the process ID, which
is also the first part of the subdirectory in which the process is run (the full
subdirectory will be something like `32/9124a1dj56n2346236245i2343`, so just a
longer hash). We then get the process and its name. Lastly, we get how many
instances of each process are currently running or have finished. Here we only
have the one process, of course, but this will soon change.

* Let's check that everything worked: type `ls results/` and see that it
  contains the output we expected.

* Let's explore the working directory: change into whatever directory is
  specified by the process ID (your equivalent to `work/32/9124a1[...]`).

What do you see when you list the contents of this directory? You should see a
symbolic link named `a.txt` pointing to the real location of this file, plus a
normal file `a.upper.txt`, which is the output of the process that was run in
this directory. You generally only move into these work directories when
debugging errors in your workflow, and Nextflow has some tricks to make this
process a lot easier - more on this later.

So, in summary: we have three components: a set of inputs stored in a channel, a
set of processes and a workflow that defines which processes should be run in
what order. We tell Nextflow to *push* the inputs through the entire workflow,
so to speak.

* Now it's your turn! Move back to the workflow root and make it use only the
  `b.txt` input file and give you the `b.upper.txt` instead.

* Run your workflow and make sure it works before you move on; check below if
  you're having trouble.

<details>
<summary> Click to show </summary>

```nextflow
ch_input = Channel.fromPath( "b.txt" )
```

</details>

# Viewing channel contents

Something that's highly useful during development of Nextflow workflows is to
view the contents of channels, which can be done with the `view()` operator.

* Add the following to your workflow definition (on a new line) and execute the
  workflow: `ch_input.view()`. What do you see?

* Remove the `view()` operator once you're done.

It can be quite helpful to view the channel contents whenever you're unsure of
what a channel contains or if you've run into some kind of bug or error, or even
just when you're adding something new to your workflow. Remember to view the
channel contents whenever you need to during the rest of this tutorial!

# Files and sample names

One powerful feature of Nextflow is that it can handle complex data structures
as input, and not only filenames. One of the more useful things this allows us
to do is to couple sample names with their respective data files inside
channels.

* Change the channel definition to the following:

```nextflow
ch_input = Channel
    .fromPath ( "a.txt" )
    .map      { file -> tuple(file.getBaseName(), file) }
```

Here we create a *tuple* (something containing multiple parts) using the `map`
operator, the *base name* of the file (`a`) and the file path (`a.txt`). The
statement `.map{ file -> tuple(file.getBaseName(), file) }` can thus be read as
"replace the channel's contents with a tuple containing the base name and the
file path". The contents of the channel thus change from `[a.txt]` to `[a,
a.txt]`. Passing the sample name or ID together with the sample data in this way
is extremely useful in a workflow context and can grealy simplify downstream
processes.

Before this will work, however, we have to change the process itself to make use
of this new information contained in the `ch_input` channel.

* Change the process definition to the following:

```nextflow
process CONVERT_TO_UPPER_CASE {
    publishDir "results/",
        mode: "copy"

    input:
    tuple val(sample), path(file)

    output:
    path("${sample}.upper.txt")

    script:
    """
    tr [a-z] [A-Z] < ${file} > ${sample}.upper.txt
    """
}
```

Notice how the input now is aware that we're passing a tuple as input, which
allows us to use both the `file` variable (as before) and the new `sample`
variable. All that's left now is to change the input to our pipeline!

* Change the channel definition line from `.fromPath ( "a.txt" )` to
  `.fromPath ( ["a.txt", "b.txt"] )` and try running the pipeline. Make sure it
  works before you move on! Remember to use the `view()` operator if you want to
  inspect the channel contents in detail.

# Input from samplesheets

So far we've been specifying inputs using strings inside the workflow itself,
but hard-coding inputs like this is not ideal. A better solution is to use
samplesheets instead, *e.g.* comma- or tab-separated data files; this is
standard for many pipelines, including [nf-core](https://nf-co.re/). Take, for
example, the following CSV file:

```no-highlight
a,a.txt
b,b.txt
```

This specifies the samples and their respective files on each row. Using such a
file is much more portable, scalable and overall easier to use than simple
hard-coding things in the workflow definition itself. We might also include an
arbitrary number of additional metadata columns, useful for downstream
processing and analyses. Using contents of files as input can be done using the
`.splitCsv()` and `.map{}` operators, like so:

```nextflow
ch_input = Channel
    .fromPath ( "first_samplesheet.csv" )
    .splitCsv ( )
    .map      { row -> tuple(row[0], file(row[1])) }
```

The `.SplitCsv()` operator lets the channel know the input is a CSV file, while
the `.map{}` operator makes the CSV content into a tuple from the first and
second elements of each row.

* Change the input channel definition to the code above and create the
  `first_samplesheet.csv` file as shown above.

* Add the `.view()` operator somewhere to show the contents of `ch_input`.

* Execute the pipeline. Do you see what you expect? Remove the `.view()`
  operator before moving on.

> **Note** <br>
> While we are still hard-coding the name of the samplesheet it is still much
> better to edit a samplesheet than having to edit the pipeline itself - there
> are also convenient ways to work around this using *parameters*, which we'll
> talk more about later in this tutorial.

We can also specify a header in our samplesheet like so: `.splitCsv(header:
true)`.

* Add an appropriate header to your samplesheet, make sure your workflow can
  read it and execute. Use `.view()` to see what's going on, if needed.

# Adding more processes

It's time to add more processes to our workflow! We have the two files
`a.upper.txt` and `b.upper.txt`; the next part of the workflow is a step
that concatenates the content of all these UPPERCASE files.

We already have a channel containing the two files we need: the output of the
`CONVERT_TO_UPPER_CASE` process called `CONVERT_TO_UPPER_CASE.out`. We can use
this output as input to a new process using the syntax:
`CONVERT_TO_UPPER_CASE.out.collect()`. The `collect()` operator groups all the
outputs in the channel into a single data object for the next process. This is a
*many-to-one* type of operation: a stream with several files (*many*) is merged
into a lone list of files (*one*). If `collect()` was not used, the next process
would try to run a task for each file in the output channel.

Let's put this in use by adding a new process to the workflow definition. We'll
call this process `CONCATENATE_FILES` and it will take the output from
`CONVERT_TO_UPPER_CASE` as input, grouped using the `collect()` operator.

* Add a line to your workflow definition for this new process with the
  appropriate input - remember that you can use `.view()` to check channel
  contents; click below if you're having trouble.

<details>
<summary> Click to show </summary>

```nextflow
CONCATENATE_FILES( CONVERT_TO_UPPER_CASE.out.collect() )
```

</details>

Now all we have to do is define the actual `CONCATENATE_FILES` process in the
process definition section.

* Copy the following code as a new process into your workflow:

```nextflow
process CONCATENATE_FILES {
    publishDir "results/",
        mode: "copy"

    input:
    path(files)

    output:
    path("*.txt")

    script:
    """
    cat ${files} > concat.txt
    """
}
```

* Run your workflow again and check the `results/` directory. At this point you
  should have three files there: `a.upper.txt`, `b.upper.txt` and `concat.txt`.

* Inspect the contents of `concat.txt` - do you see everything as you expected?

Note the use of `path(files)` as input. Although we pass a list of files as
input, the list is considered a single object, and so the `files` variable
references a list. Each file in that list can be individually accessed using an
index e.g. `${files[0]}`, or as we do here, use the variable without an index
to list all the input files.

> **Quick recap** <br>
> In this section we've learnt:
>
> * How to create, execute and extend workflows
> * How to explore the `work` directory and channel contents
> * How to couple sample names to sample data files
> * How to use samplesheets as input
> * How to collect multiple files as single inputs for processes
