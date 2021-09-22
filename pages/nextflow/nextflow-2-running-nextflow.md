Before we go into the details regarding how to write a workflow in Nextflow,
let's just try running the MRSA workflow! The simplest and most common way of
running a Nextflow workflow is by `nextflow run main.nf`, which means "run the
workflow situated in the current directory". Make sure you've activated the
`nextflow-env` Conda environment and that you are standing in `workshop-reproducible-research/tutorials/nextflow`
and run the command. After it has finished (which shouldn't take long), you
should see something like this:

```bash
N E X T F L O W  ~  version 21.04.0
Launching `./main.nf` [friendly_bhaskara] - revision: b4490b9201
executor >  local (17)
[c9/e5f818] process > get_sra_by_accession (SRR935092) [100%] 3 of 3 ✔
[d5/b5f24e] process > run_fastqc (SRR935092)           [100%] 3 of 3 ✔
[91/2cea54] process > run_multiqc                      [100%] 1 of 1 ✔
[e0/b4fd37] process > get_genome_fasta                 [100%] 1 of 1 ✔
[87/32ce10] process > index_genome                     [100%] 1 of 1 ✔
[56/e9a460] process > align_to_genome (SRR935092)      [100%] 3 of 3 ✔
[ed/d8c223] process > sort_bam (SRR935092)             [100%] 3 of 3 ✔
[e7/4a6bda] process > get_genome_gff3                  [100%] 1 of 1 ✔
[e9/84f093] process > generate_counts_table            [100%] 1 of 1 ✔
```

The first few lines are information about this particular run, including the
Nextflow version used, which workflow definition file was used (`main.nf`), a
randomly generated run name (an adjective and a scientist), the revision number
as well as the executor used (local, in this case).

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

You can find the results of the workflow in the `results/` directory, while the
hidden directories where each process is run is found in `work/`. If you want to
examine a process more carefully, simply navigate to its directory, *e.g.*
`work/` followed by the first part given in the list above, *e.g.* `87/32ce10`.
If you want to view the logs you can use `nextflow log <run name>` or look in
the hidden file `.nextflow.log`.

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
> - How to execute Nextflow workflows
> - How to re-run only changed parts of workflows
> - How to get execution reports, timelines and DAG visualisations
