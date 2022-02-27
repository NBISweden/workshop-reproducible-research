We've so far been working with a relatively non-generalised work: it's got
hard-coded inputs, paths and genome references. This is perfectly fine for a
project that is purely aimed at getting reproducible results (which is the full
extent of what you want in a lot of cases), but it can be made a lot more
generalisable. Let's go through the MRSA workflow and see what can be improved!

# Configuration basics

One of the things that allow generalisability of Nextflow workflows is
*parameters*, which hold information and values that can be changed directly on
the command-line at the time of execution. One use of parameters in our MRSA
workflow is to remove the hard-coded `results` output directory, for example.
Parameters can be written on the following form:

```nextflow
params {
    parameter_1 = "some/data/path"
    parameter_2 = 42
    parameter_3 = ["a", "list", "of", "strings"]
}
```

You would then refer to these parameters using *e.g.* `params.parameter_1`
anywhere you need to in the workflow. The parameters are not put into the
`main.nf` file, but rather into a *configuration file*. The default name of this
file is `nextflow.config`.

* Create a configuration file and add a parameter for the `results` output
  directory.

* Use your newly created parameter in the `publishDir` directory of a process
  (it'll be on the form of `${params.resultsdir}/some/other/path`, for example).
  Run your workflow to see if it worked.

> **Tip** <br>
> Instead of manually changing all of the hard-coded directories in your
> workflow you can use the following little `sed` command, which will do it for
> you: `sed 's/\"results\//\"${params.resultsdir}\//g' main.nf > tmp; mv tmp
> main.nf`.

# Command line parameters

Parameters can be changed on-the-fly when executing workflows like so: `nextflow
run main.nf --parameter_name some_value`

* Run your workflow using the parameter you previously created, but pick
  something other than the default value!

You should now have a new directory containing all the results! This is highly
useful if you want to keep track of separate runs of a workflow with different
software parameters, for example: `nextflow run main.nf --important_param value1
--resultsdir value1`, or simply want to keep the results of separate versions of
the same workflow. You can, of course, also change parameters directly in the
configuration file, rather than on the command line!

* Add another parameter for the input SRA IDs and execute your workflow to check
  that it worked.

# Other configuration scopes

There are lots of things that you might want to add to you configuration, not
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

> **Quick recap** <br>
> In this section we learnt:
>
> * How to create parameters in a configuration file
> * How to specify parameters on the command line
> * How to add  workflow manifest and other configuration scopes
