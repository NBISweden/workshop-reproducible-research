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
// Enable DSL2 functionality
nextflow.enable.dsl = 2

// Workflow definition
workflow {
    // Define input files
    ch_input = Channel.fromPath("a.txt")

    // Run workflow
    CONVERT_TO_UPPER_CASE(ch_input)
}

// Process definition
process CONVERT_TO_UPPER_CASE {
    publishDir "results/", mode: "copy"

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

Here we have three separate parts. The first part enables the DSL2 (*Domain
Specific Language 2*) functionality, and is required to use some of the newer
and more powerful features of Nextflow. The next part is the *workflow
definition*, while the last is a *process*. Let's go through the last two in
more detail!

> **Nextflow comments** <br>
> Double-slashes (`//`) are used for comments in Nextflow.

> **Nextflow and whitespace** <br>
> Nextflow is not indentation-sensitive. In fact, Nextflow doesn't care at all
> about whitespace, so go ahead and use it in whatever manner you think is
> easiest to read and work with! Do keep in mind that indentations and other
> types of whitespace *does* improve readability, so it's generally not a good
> idea to forego it entirely, even though you can.

# Workflow definitions

The workflow definition here has two parts, each doing an important job for any
Nextflow workflow. The first part defines a *channel*, which is an asynchronous
first-in-first-out stream of data that connect a workflow's various inputs and
outputs. In this particular case, we define a `Channel` using the `.fromPath`
channel factory on the specific file path `a.txt`, and name the channel `ch_input`. You
can read this as *"create the channel `ch_input` and send the file
`a.txt` into it"*.

> **Naming channels** <br>
> A channel can be named anything you like, but it is good practice to prepend
> them with `ch_`, as that makes it clear which variables are channels and which
> are just normal variables.

How do we use these channels then? Channels pass data to and from processes
through our workflow. By providing channels as arguments to processes, we describe
how we want data to flow. This is exactly what we do in the second part: we call
our `CONVERT_TO_UPPER_CASE` process with the `ch_input` as input argument - this
is very similar to functional programming.

This is our entire workflow, for now: the creation of a channel followed by
using the contents of that channel as input to a single process. Let's look at
how processes themselves are defined!

# Process definitions

Looking at the process in the code above, we can see several parts. The process
block starts with its name, in this case `CONVERT_TO_UPPER_CASE`, followed by
several sections: `publishDir`, `input`, `output` and `script`.

> **Naming processes** <br>
> A process can be named using any case, but a commonly used convention is to use
> UPPERCASE letters for processes to visually distinguish them in the workflow.
> You do not have to follow this if you don't want to, but we do so here.

Let's ignore the first section for now and focus on the last three. The `input` and
`output` sections describe the data expected to come through the channel for this
specific process. Each line of `input` describes the data expected for each process
argument, in the order used in the workflow. In this case, `CONVERT_TO_UPPER_CASE`
expects a single channel (one line of input), and expects the data to be filenames
(of type `path`). Notice that there is a difference between how the inputs and
outputs are declared? The `output` is an explicit string (*i.e* surrounded by
quotes), while the input is a variable named `file`. This means inputs can be
referenced in the process without naming the data explicitly, unlike the output
where the name needs to be explicit. We'll get back to exactly how
this works in just a moment.

Let's move on to the first section: `publishDir`. This tells
Nextflow where the output of the process should be stored when it is finished;
setting `mode` to `"copy"` just means that we want to copy the output files to
the publishing directory, rather than using a symbolic link (which is the
default).

# Executing workflows

Let's try running the workflow we just created!

* Type the following in your terminal:

```bash
nextflow run main.nf
```

This will make Nextflow run the workflow specified in your `main.nf` file. You
should see something along these lines:

```no-highlight
N E X T F L O W  ~  version 21.04.0
Launching `./main.nf` [mad_legentil] - revision: 87f0c253ed
executor >  local (1)
[32/9124a1] process > CONVERT_TO_UPPER_CASE (1) [100%] 1 of 1 ✔
```

The first few lines are information about this particular run, including the
Nextflow version used, which workflow definition file was used, a randomly
generated run name (an adjective and a scientist), the revision ID as well
as where the processes were executed (locally, in this case).

What follows next is a list of all the various processes for this particular
workflow. The order does not necessarily reflect the order of execution
(depending on each process’ input and output dependencies), but they are in the
order they were defined in the workflow file - there's only the one process
here, of course. The first part (*e.g* `[32/9124a1]`) is the process ID, which
is also the first part of the subdirectory in which the process is run (before
the outputs are transferred to the publish directory). We then get the process
and its name. Lastly, we get how many instances of each process are currently
running or have finished. Here we only have the one process, of course, but this
will soon change.

* Let's check that everything worked: type `ls results/` and see that it
  contains the output we expected.

* Let's explore the working directory: change into whatever directory is
  specified by the process ID (your equivalent to `work/32/9124a1[...]`).

What do you see when you list the contents of this directory? You should,
hopefully, see a symbolic link named `a.txt` pointing to the real location of
this file, plus a normal file `a.upper.txt`, which is the output of the process
that was run in this directory. While it seems cumbersome to manually move into
these work directories it is something you only do when debugging errors in your
workflow, and Nextflow has some tricks to make this process a lot easier - more
on this later.

So, how does this all work? Well, we have three components: a set of inputs, a
set of processes and a workflow that defines which processes should be run. We
tell Nextflow to *push* the inputs through the entire workflow, so to speak.

* Now it's your turn! Move back to the workflow root and make it use only the
  `b.txt` input file and give you the `b.upper.txt` instead.

* Run your workflow and make sure it works before you move on.

# Files and sample names

Having to manually change inputs and outputs like you just did is not really
ideal, is it? Hard-coding outputs is rarely good, so let's try to change that.
One powerful feature of Nextflow is that it can handle complex data structures
as input, and not only filenames. One strategy we can follow is to create
a prefix for our output and pass it together with the filename.

* Change the channel definition to the following:

```nextflow
ch_input = Channel
    .fromPath("a.txt")
    .map{ file -> tuple(file.getBaseName(), file) }
```

Okay, so what does that do, exactly? Well, the added line containing the
`.map{}` statement changes the data stream to be `[prefix, file]` instead of
just `[file]` - we generate the prefix from the *base name* of the file itself,
*i.e.* the file without extension or directory. We now have to change the
process itself to make use of this new information contained in the `ch_input`
channel.

* Change the process definition to the following:

```nextflow
process CONVERT_TO_UPPER_CASE {
    publishDir "results/", mode: "copy"

    input:
    tuple val(prefix), path(file)

    output:
    path("${prefix}.upper.txt")

    script:
    """
    tr [a-z] [A-Z] < ${file} > ${prefix}.upper.txt
    """
}
```

Notice how the input now is aware that we're passing a *tuple* as input, which
allows us to use both the `file` variable (as before) and the new `prefix`
variable. All that's left now is to change the input to our pipeline!

* Change the channel definition line from `.fromPath("a.txt")` to
  `.fromPath(["a.txt", "b.txt"])` and try running the pipeline. Make sure it
  works before you move on!

# Adding more processes

It's time to add more processes to our workflow! We have the two files
`a.upper.txt` and `b.upper.txt`; the next part of the workflow is a step
that concatenates the content of all these UPPERCASE files.

We already have a channel containing the two files we need: the output of the
`CONVERT_TO_UPPER_CASE` process called `CONVERT_TO_UPPER_CASE.out`. We can use
this output as input to a new process using the syntax:
`CONVERT_TO_UPPER_CASE.out.collect()`. The `collect()` operator, groups all the
outputs in the channel into a single data object for the next process. This is
a *many-to-one* type of operation: a stream with several files (*many*) is
merged into a lone list of files (*one*). If `collect()` was not used, the next
process would try to run a task for each file in the output channel.

Let's put this in use by adding a new process to the workflow definition. We'll
call this process `CONCATENATE_FILES` and it will take the output from
`CONVERT_TO_UPPER_CASE` as input, grouped using the `collect()` operator.

* Add a line to your workflow definition for this new process with the
  appropriate input - click below if you're having trouble.

??? example "Click to show the solution"
    ```nextflow
    CONCATENATE_FILES( CONVERT_TO_UPPER_CASE.out.collect() )
    ```

Now all we have to do is define the actual `CONCATENATE_FILES` process in the
process definition section.

* Copy the following code as a new process into your workflow:

```nextflow
process CONCATENATE_FILES {
    publishDir "results/", mode: "copy"

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

Run your workflow again and check the `results/` directory. At this point you
should have three files there: `a.upper.txt`, `b.upper.txt` and `concat.txt`.

* Inspect the contents of `concat.txt` - do you see everything as you expected?

Note the use of `path(files)` as input. Although we pass a list of files as
input, the list is considered a single object, and so the `files` variable
references a list. Each file in that list can be individually accessed using an
index e.g. `${files[0]}`, or as we do here, use the variable without an index
to list all the input files.

# Viewing channel contents

As our channels become more complicated it is useful to actually check out
what's inside them: you can do this using the `.view()` operator.

* Add the following to your workflow definition (on a new line) and execute the
  workflow: `ch_input.view()`. What do you see?

It can be quite useful to inspect channel contents like this when you are
developing workflows, especially if you are working with tuples, maps and any
transforming operators in general.

* Check the channel contents of the (1) raw and (2) collected output of the
  `CONVERT_TO_UPPER_CASE` process. How are they different?

!!! Success "Quick recap"
    In this section we've learnt:

    * How to create and extend simple Nextflow workflows
    * How to create channels for input data
    * How to execute workflows
    * How to explore Nextflow's `work` directory
    * How to generalize workflows
    * How to view channel contents
