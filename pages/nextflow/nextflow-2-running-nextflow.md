Before we go into the details regarding how to write a workflow in Nextflow,
let's just try running the MRSA workflow! The simplest and most common way of
running a Nextflow workflow is by `nextflow run main.nf`, which means "run the
workflow defined in the `main.nf` file".

* Make sure you've activated the `nextflow-env` Conda environment and that you
  are standing in `workshop-reproducible-research/tutorials/nextflow` and run
  the `nextflow run main.nf` command.

After it has finished (which shouldn't take long), you should see something
like this:

```bash
N E X T F L O W  ~  version 21.04.0
Launching `./main.nf` [friendly_bhaskara] - revision: b4490b9201
executor >  local (17)
[c9/e5f818] process > GET_SRA_BY_ACCESSION (SRR935092) [100%] 3 of 3 ✔
[d5/b5f24e] process > RUN_FASTQC (SRR935092)           [100%] 3 of 3 ✔
[91/2cea54] process > RUN_MULTIQC                      [100%] 1 of 1 ✔
[e0/b4fd37] process > GET_GENOME_FASTA                 [100%] 1 of 1 ✔
[87/32ce10] process > INDEX_GENOME                     [100%] 1 of 1 ✔
[56/e9a460] process > ALIGN_TO_GENOME (SRR935092)      [100%] 3 of 3 ✔
[ed/d8c223] process > SORT_BAM (SRR935092)             [100%] 3 of 3 ✔
[e7/4a6bda] process > GET_GENOME_GFF3                  [100%] 1 of 1 ✔
[e9/84f093] process > GENERATE_COUNTS_TABLE            [100%] 1 of 1 ✔
```

The first few lines are information about this particular run, including the
Nextflow version used, which workflow definition file was used, a randomly
generated run name (an adjective and a scientist), the revision number as well
as the executor used (local, in this case).

What follows next is a list of all the various processes for this particular
workflow. Just like Snakemake, the order does not necessarily reflect the order
of execution (depending on each process' input and output dependencies), but
they are in the order they were defined in the workflow file. The first part
(*e.g* `[c9/e5f818]`) is the process ID, which is also the first part of the
subdirectory in which the process is run. We then get the process and its name,
where some include the current sample name being run (*e.g.* `SRR935092`).
Lastly, we get how many instances of each process are currently being and have
been run.

If you successfully ran the command above, you will have seen that the output
shown above is not exactly what it looks like during execution: the list of
processes is continuously updated and changed, depending on what processes are
being run and which ones have been completed. You might also get a different SRA
ID than above as the "final" one, *i.e.* the last one being handled by the
different processes. This is because Nextflow (just like Snakemake) doesn't
order the execution of the processes by default, so it happens randomly.

You can find the results of the workflow in the `results/` directory, with a
directory structure specified in the workflow itself (we'll come to that in the
next part of this tutorial).

* Check what's inside the results directory by running `ls results/`.

The directories you see should be familiar to you if you went through the
Snakemake tutorial, as they are the same.

You can find the work done by each process in the `work/` directory. If you want
to examine a process more carefully, simply navigate to its directory, *e.g.*
`work/` followed by the first part given in the list above, *e.g.* `87/32ce10`.

* Let's check the results of the `RUN_MULTIQC` process! Move to the directory
  corresponding to the `work/91/2cea54(...)` output in your computer and list
  the files using `ls -l`.

You should see a number of things: the `multiqc.html` report, a `multiqc_data`
directory and three symbolic links to files ending in `.zip`. Why are some files
symbolic links, rather than normal files? This is because Nextflow runs each
process in its own directory, using symbolic links for files coming from
previously run processes. Think about the `RUN_MULTIQC` process: it takes the
output from the `RUN_FASTQ` process (*i.e.* the compressed files) and runs the
MultiQC software on them. The HTML report and the data directory are outputs of
the `RUN_MULTIQC` process, so they are created here. If there was any downstream
process that would use any of these files, they would be symbolically linked
from here!

> **Nextflow logs** <br>
> If you want to view the logs you can use `nextflow log <run name>` or look in
> the hidden file `.nextflow.log`. They are quite detailed and can greatly
> facilitate debugging your workflow when you run into errors!

## Re-running workflows

Just like in Snakemake, something you often want to do in Nextflow is to re-run
the workflow when you changed some input files or some of the code for its
analyses, but you don't want to re-run the entire workflow from start to finish.
Let's find out how this works in Nextflow!

* Run the same `nextflow run main.nf` command again.

What happened here? Nextflow actually re-ran the entire workflow from scratch,
even though we didn't change anything. This is the default behaviour of
Nextflow.

* Let's try that again: run `nextflow run main.nf -resume` instead.

Now you can see that Nextflow didn't actually re-run anything, but used the
cached results from the previous run by supplying the `-resume` flag to the
command!

As mentioned in the introduction, how Nextflow re-runs its processes is
slightly different than Snakemake, in that Snakemake only checks if the inputs
are newer than the outputs, whereas Nextflow also checks if any of the code or
scripts have changed. This means that when you change the `script` section of a
Snakemake rule (or an external script it calls) you need to run with the `-R
<rule>` flag in order to re-run with the new changes, whereas Nextflow
automatically keeps track of code changes for you.

## Reports and visualisations

Nextflow has a number of reports and visualisations it can provide you
automatically for any given pipeline you have. Let's have Nextflow create a few
of them for us by executing the following command:

```bash
nextflow run main.nf -with-report -with-timeline -with-dag dag.png
```

After successful executing, you will find three more files in your current
directory: `report.html`, `timeline.html` and `dag.png`. The first file contains
a workflow report, which includes various information regarding execution such
as runtime, resource usage and details about the different processes. The second
file contains a timeline for how long each individual process took to execute,
while the last contains a visualisation of the workflow itself.

Take a few minutes to browse these files for yourself! When running a workflow
you can of course choose which of these additional files you want to include by
picking which ones are important or interesting to you - or don't include any!

> **Quick recap:** <br>
> In this section we covered:
>
> - How to execute and examine the output of Nextflow workflows
> - How to re-run only changed parts of workflows
> - How to get execution reports, timelines and DAG visualisations
