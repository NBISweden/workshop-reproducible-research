It's time to start working with a more realistic workflow using the MRSA case
study of this course! We've created a bare-bones version of this pipeline,
but we'll work our way through it as we go along and learn more about Nextflow's
features and functionality. We'll begin with some specifics of how workflows are
executed and what you can get from them.

# Reports and visualisations

Let's start with running the workflow plus getting some reports and
visualisation while we're at it!

* Run the workflow using the following command: `nextflow run main_mrsa.nf
  -with-report -with-timeline -with-dag dag.png`.

After successful executing, you will find three more files in your current
directory: `report.html`, `timeline.html` and `dag.png`. The first file contains
a workflow report, which includes various information regarding execution such
as runtime, resource usage and details about the different processes. The second
file contains a timeline for how long each individual process took to execute,
while the last contains a visualisation of the workflow itself.

Take a few minutes to browse these files for yourself! When running a workflow
you can of course choose which of these additional files you want to include by
picking which ones are important or interesting to you - or don’t include any!

# Re-running workflows

Just like in Snakemake, something you often want to do in Nextflow is to re-run
the workflow when you changed some input files or some of the code for its
analyses, but you don't want to re-run the entire workflow from start to finish.
Let’s find out how this works in Nextflow!

* Run the same `nextflow run main_mrsa.nf` command again.

What happened here? Nextflow actually re-ran the entire workflow from scratch,
even though we didn't change anything. This is the default behaviour of
Nextflow.

* Let’s try that again: `nextflow run main_mrsa.nf -resume` instead.

Now you can see that Nextflow didn't actually re-run anything. The `-resume` 
flag instructed Nextflow to use the cached results from the previous run!

As mentioned in the introduction, how Nextflow re-runs its processes is slightly
different from Snakemake, in that Snakemake only checks if the inputs are newer
than the outputs, whereas Nextflow also checks if any of the code or scripts
have changed. This means that when you change the script section of a Snakemake
rule (or an external script it calls) you need to run with the `-R <rule>` flag
in order to re-run with the new changes, whereas Nextflow automatically keeps
track of code changes for you.

* Remove the `results` directory and re-run the workflow again using the
  `-resume` flag.

We removed all the results we used before, but we still managed to resume the
workflow and use its cache - how come? Remember that Nextflow uses the `work`
directory to run all of its tasks, while the `results` directory is just where
we have chosen to publish our outputs. We can thus delete the `results`
directory all we want and still get everything back without having to re-run
anything. If we were to delete the `work` directory, however...

* Delete the `work` directory and re-run the workflow using the `-resume` flag.

There is no longer any cache for Nextflow to use, so it re-runs from the start!
This is good to keep in mind: you can always delete the output directories of
your workflow, but if you mess with `work` you'll lose, well... work!

# Logs

Nextflow keeps a log of all the workflows that have been executed. Let's check
it out!

* Type `nextflow log` to get a list of all the executions.

Here we get information about when the workflow was executed, how long it ran,
its run name, whether it succeeded or not and what command was used to run it.
You can also use `nextflow log <run name>` to show the work subfolder of each 
task that was executed for that run. We can, however, get even more detailed 
information about the latest run by looking into the `.nextflow.log` file!

* Look into the latest log by typing `less .nextflow.log`.

You'll be greeted by a wealth of debugging information, which may even seem a
bit overkill at this point! This level of detail is, however, quite useful both
as a history of what you've attempted and as an additional help when you run 
into errors! Also, it helps with advanced debugging - which we'll get into later.

> **Quick recap** <br>
> In this section we've learnt:
>
> * How to get automatic reports and visualisations
> * How to re-run workflows
> * How to check the Nextflow logs
