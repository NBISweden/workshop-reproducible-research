We'll start by creating a very simple workflow from scratch, to show how
Nextflow works. The workflow will be the same as in the Snakemake tutorial: it
will take two input files and convert them to UPPERCASE letters.

* Start by running the following commands:

```bash
touch main.nf
echo "This is a.txt" > a.txt
echo "This is b.txt" > b.txt
```

Open the `main.nf` file with an editor of your choice. This is the main workflow
file used in Nextflow, analogous to the `Snakefile` in Snakemake.

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

This looks a bit more complicated than the single rule we had in the Snakemake
tutorial for the same workflow, but there are reasons for that that we'll
explain in a moment. Here we have three separate parts. The first part enables
the DSL2 (*Domain Specific Language 2*) functionality, and is required to use
some of the newer and more powerful features of Nextflow. The next part is the
workflow definition, while the last is a *process*, equivalent to Snakemake's
rules. Let's go through the last two in more detail!

> **Nextflow comments** <br>
> Double-slashes (`//`) are used for comments in Nextflow.

> **Nextflow and whitespace** <br>
> Unlike Snakemake, Nextflow is not indentation-sensitive. In fact, Nextflow
> doesn't care at all about whitespace, so go ahead and use it in whatever
> manner you think is easiest to read and work with! Do keep in mind that
> indentations and other types of whitespace *does* improve readability, so it's
> generally not a good idea to forego it entirely, even though you can.

## Process definitions

There are quite a few similarities between Snakemake's rules and Nextflow's
processes, which can be seen by just a cursory glance. We have both `input`,
`output` and `script` (`shell` in Snakemake) parts, for example, and we
similarly have named processes (here `CONVERT_TO_UPPER_CASE`). The `script` part
is what will be executed when the process is run, exactly as in Snakemake, but
the `input` and `output` parts are slightly different.

> **Naming processes** <br>
> A process can be named anything you like, but there exists a commonly used
> convention to use UPPERCASE letters for processes. You do not have to follow
> this if you don't want, but we do so here.

Notice that there is a difference between how the inputs and outputs are
declared? The `output` is an explicit string (*i.e* surrounded by quotes), while
the input is a variable named `file`. The reason for this is that processes are
like functions, meaning they can have have varying inputs (and arguments in
general). There is one part of the process definition that has no equivalent in
Snakemake, and that is the `publishDir` directive. What this does is tell
Nextflow where the output of the process should be stored when it is finished;
setting `mode` to `"copy"` just means that we want to copy the output files to
the publishing directory, rather than using a symbolic link (which is the
default).

Now that we know the basics of processes, let's move on to the workflow
definition!

## Workflow definitions

The workflow definition here has two parts, each doing an important job for any
Nextflow workflow. The first part defines a *channel*, which is an asynchronous
first-in-first-out stream of data that connect a workflow's various inputs and
outputs. In this particular case, we define a `Channel` using the `.fromPath`
operator on the specific file path `a.txt` and name the channel `ch_input`. This
means that `ch_input` contains a stream with the file `a.txt`.

> **Naming channels** <br>
> A channel can be named anything you like, but it is good practice to prepend
> them with `ch_`, as that makes it clear which variables are channels and which
> are just normal variables.

What can we do with these channels, then? The idea is that we put the channel
through our workflow, *i.e.* the various processes we want to run on the inputs
stored in the channel. This is exactly what we do in the second part: we call
our `CONVERT_TO_UPPER_CASE` process with the `ch_input` as input argument.

## Executing workflows

Let's try running the workflow we just created!

* Type the following in your terminal:

```bash
nextflow run main.nf
```

This will make Nextflow run the workflow specified in your `main.nf` file. You
should see something along these lines:

```no-highlight
N E X T F L O W  ~  version 21.10.4
Launching `./main.nf` [mad_legentil] - revision: 87f0c253ed
executor >  local (1)
[32/9124a1] process > CONVERT_TO_UPPER_CASE (1) [100%] 1 of 1 ✔
```

The first few lines are information about this particular run, including the
Nextflow version used, which workflow definition file was used, a randomly
generated run name (an adjective and a scientist), the revision number as well
as how the run was executed (locally, in this case).

What follows next is a list of all the various processes for this particular
workflow. Just like Snakemake, the order does not necessarily reflect the order
of execution (depending on each process’ input and output dependencies), but
they are in the order they were defined in the workflow file - there's only the
one process here, of course. The first part (*e.g* `[32/9124a1]`) is the process
ID, which is also the first part of the subdirectory in which the process is
run. We then get the process and its name. Lastly, we get how many instances of
each process are currently being and have been run. Here we only have the one
process, of course, but this will soon change.

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
This is the opposite of Snakemake, where we define a workflow and then ask
Snakemake to give us a final output we know it can create - *pull*.

* Now it's your turn! Move back to the workflow root and change the it to use
  the `b.txt` input file and give you the `b.upper.txt` instead.

* Run your workflow and make sure it works before you move on.

## Files and sample names

Having to manually change inputs and outputs like you just did is not really
ideal, is it? Hard-coding things is rarely good, so let's try to change that.
What we need is a way to use the sample name for each input file, something that
does something similar to *wildcards* in Snakemake. Here is where one of the
features of Nextflow really shines: being able to parse values *alongside* files
in channels!

* Change the channel definition to the following:

```nextflow
ch_input = Channel
    .fromPath("a.txt")
    .map{file -> tuple(file.getBaseName(), file)}
```

Okay, so what does that do, exactly? Well, the added line containing the
`.map{}` statement changes the input stream to be `[name, file]` instead of just
`[file]` - we get the name from the *base name* of the file itself, *i.e.* the
file without extension or directory. We now have to change the process itself to
make use of this new information contained in the `ch_input` channel.

* Change the process definition to the following:

```nextflow
process CONVERT_TO_UPPER_CASE {
    publishDir "results/", mode: "copy"

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

Notice how the input now is aware that we're passing a *tuple* as input, which
allows us to use both the `file` variable (as before) and the new `sample`
variable. All that's left now is to change the input to our pipeline!

* Change the channel definition line from `.fromPath("a.txt")` to
  `.fromPath(["a.txt", "b.txt"])` and try running the pipeline. Make sure it
  works before you move on!

## Adding more processes

It's time to add more processes to our workflow! Here is where the differences
between Snakemake and Nextflow become more apparent: at this stage in our
Snakemake workflow we'd have two files `a.upper.txt` and `b.upper.txt` stored in
explicit directories (as specified by the pipeline). This is not the case for
Nextflow: here we have a *channel* that contains those same files, but we don't
have to care about their location. The next part of the workflow, if you recall,
is a step that concatenates the content of all the UPPERCASE files.

We already have a channel containing the two files we need: the output of the
`CONVERT_TO_UPPER_CASE` process. We can use this output as input to a new
process (which we'll name `CONCATENATE_FILES`) like so: `CONCATENATE_FILES(
CONVERT_TO_UPPER_CASE.out.collect() )`. There are two new things here, the first
of which is the `.out` attribute. As you might imagine, this refers to the
output of a particular channel. The `collect()` operator, on the other hand,
*collects* all the outputs into a single entry. This is a *many-to-one* type of
operation: a stream with several files (*many*) is merged into a lone list of
files (*one*).

As our channels become more complicated it is useful to actually check out
what's inside them: you can do this using the `view()` operator.

* Add the following to your workflow definition (on a new line) and execute the
  workflow: `ch_input.view()`. What do you see?

* Do the same for the (1) raw and (2) collected output of the `CONVERT_TO_UPPER`
  process. How are they different?

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
    cat ${files} > a_b.txt
    """
}
```

Now all we have to do is to add this process to the workflow definition, so that
it will actually be run during execution!

* Add a line to your workflow definition for this new process with the
  appropriate input - click below if you're having trouble.

<details>
<summary> Click to show </summary>

```nextflow
CONCATENATE_FILES( CONVERT_TO_UPPER_CASE.out.collect() )
```

</details>

Run your workflow again and check the `results/` directory. At this point you
should have three symbolic links there: `a.upper.txt`, `b.upper.txt` and
`a_b.txt`.

> **Quick recap** <br>
> In this section we've learnt:
>
> * How to create and extend simple Nextflow workflows
> * How to create channels for input data
> * How to execute workflows
> * How to explore Nextflow's `work` directory
> * How to generalize workflows
