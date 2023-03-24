We've so far been working with a relatively non-generalised workflow: it's got
hard-coded inputs, paths and genome references. This is perfectly fine for a
project that is purely aimed at getting reproducible results (which is the full
extent of what you want in a lot of cases), but it can be made a lot more
generalisable. Let's go through the MRSA workflow and see what can be improved!

# Parameters

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
`main_mrsa.nf`, it is preferable to define them in a separate *configuration
file*. The default name of this file is `nextflow.config` and if such a file is
present it will be used automatically by Nextflow (to supply a config file with
another name use `nextflow -c <path-to-config-file> run main_mrsa.nf`)

* Create a configuration file and add a parameter for the `results` output
  directory.

* Use your newly created parameter in the `publishDir` directory of a process
  Run your workflow to see if it worked; click below if you need help.

<details>
<summary> Click to show </summary>

```nextflow
// Configuration file
params {
    outdir = "results"
}

// A publishDir directive in a process
publishDir: "{params.outdir}", mode: "copy"
```

</details>

# Command line parameters

Workflow parameters can be assigned on the command-line by executing workflows
like so: `nextflow run main_mrsa.nf --parameter_name 'some_value'`. The workflow
parameter `parameter_name` is prefixed by a double dash `--` to tell Nextflow
this is a parameter to the workflow (a single dash is a parameter to Nextflow,
*e.g.* `-resume`). The value is also quoted (this is important for parameters
that take file paths as values).

* Run your workflow using the parameter you previously created, but pick
  something other than the default value!

You should now have a new directory containing all the results! This is highly
useful if you want to keep track of separate runs of a workflow with different
software parameters, for example: `nextflow run main.nf --important_param
'value1' --resultsdir 'results-value1'`, or simply want to keep the results of
separate versions of the same workflow. You can also change parameters by using
the `-params-file` option or by using another configuration file (and using
`-c`), rather than on the command line!

# Configuring inputs

Remember the input for the MRSA workflow, the `ch_input` channel? This input
(the `samplesheet.csv` file) is hard-coded inside the `main_mrsa.nf` file. This
could also be made into a parameter!

* Change the definition of the `ch_input` channel to take the value of a new
  parameter of your choice, defined in the configuration file.

You should now have a more generalised input to your workflow! Try to run it to
make sure it works - look below if you need some help.

<details>
<summary> Click to show </summary>

```nextflow
// Channel definition
ch_input = Channel
    .fromPath ( params.input )
    .splitCsv ( header: true )

// Configuration file
input = "samplesheet.csv"
```

</details>

By specifying inputs from sample sheets like this we can change inputs of a
workflow execution by creating another sample sheet and specifying *e.g.*,
`--input samplesheet-2.csv` on the command line. This is highly useful when you
want to run a single sample *e.g.*, when testing a workflow, or when you want to
keep track of all the different inputs you've used historically.

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

> **Quick recap** <br>
> In this section we learnt:
>
> * How to create parameters in a configuration file
> * How to specify parameters on the command line
> * How to add workflow manifest and other configuration scopes
