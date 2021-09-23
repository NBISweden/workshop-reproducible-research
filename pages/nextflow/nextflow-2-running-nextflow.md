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

If you run the same command again, Nextflow will re-run the entire workflow
from scratch. This is an important difference to Snakemake, which only re-runs
those parts of the workflow that have changed. This behaviour exists in
Nextflow as well through the `-resume` flag. Nextflow does, however, not only
keep track of changed inputs and outputs, but also processes and parameters.
This means you don't have to run with a specific flag and target if you changed
a process (like `-R <rule>` in Snakemake), but just re-run with `-resume` and
Nextflow will take care of the rest.

## Reports and visualisations

Nextflow can also automatically supply execution reports by running with the
`-with-report` flag. This will give you the `report.html` file, containing
various kinds of information regarding the execution, such as runtime,
resource usage, and details about the different processes. Similarly, you can 
also get a timeline report by using `-with-timeline` or a visualisation of the
entire DAG using `-with-dag <filename>.png`.

> **Quick recap:** <br>
> In this section we covered:
>
> - How to execute and examine the output of Nextflow workflows
> - How to re-run only changed parts of workflows
> - How to get execution reports, timelines and DAG visualisations
