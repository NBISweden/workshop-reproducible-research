If you want to read more about Snakemake in general you can find several
resources here:

* The Snakemake documentation is available on [ReadTheDocs](
  https://snakemake.readthedocs.io/en/stable/#).
* Here is another (quite in-depth) [tutorial](
  https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html#tutorial).
* If you have questions, check out [stack overflow](
  https://stackoverflow.com/questions/tagged/snakemake).

# Using containers in Snakemake

Snakemake also supports defining a Singularity or Docker container for each
rule (you will have time to work on the [Containers tutorial](containers-1-introduction)
later during the course). Analogous to using a rule-specific Conda environment,
specify `container: "docker://some-account/rule-specific-image"` in the rule
definition. Instead of a link to a container image, it is also possible to
provide the path to a `*.sif` file (= a Singularity file). When executing
Snakemake, add the `--use-singularity` flag to the command line. For the given
rule, a Singularity container will then be created from the image or Singularity
file that is provided in the rule definition on the fly by Snakemake and the
rule will be run in this container.

You can find pre-made Singularity or Docker images for many tools on
[https://biocontainers.pro/](https://biocontainers.pro/) (bioinformatics-specific)
or on [https://hub.docker.com/](https://hub.docker.com/).

Here is an example for a rule and its execution:

```python
rule align_to_genome:
    output:
        temp("results/bam/{sample_id,\w+}.bam")
    input:
        fastq = "data/{sample_id}.fastq.gz",
        index = expand("results/bowtie2/{genome_id}.{substr}.bt2",
            genome_id=config["genome_id"],
            substr=["1", "2", "3", "4", "rev.1", "rev.2"])
    log:
        expand("results/logs/align_to_genome/{{sample_id}}_{genome_id}.log",
            genome_id = config["genome_id"])
    container: "docker://quay.io/biocontainers/bowtie2:2.5.0--py310h8d7afc0_0"
    shell:
        """
        bowtie2 -x results/bowtie2/{config[genome_id]} -U {input.fastq} > {output} 2>{log}
        """
```

Start your Snakemake workflow with the following command:

```bash
snakemake --use-singularity
```

Feel free to modify the MRSA workflow according to this example. As Singularity
is a container software that was developed for HPC clusters, and for example the
Mac version is still a beta version, it might not work to run your updated
Snakemake workflow with Singularity locally on your computer.
In the next section we explain how you can run Snakemake workflows on UPPMAX
where Singularity is pre-installed.

# Running Snakemake workflows on UPPMAX

There are several options to execute Snakemake workflows on UPPMAX (a HPC
cluster with the SLURM workload manager). In any case, we highly recommend to
use a session manager like [tmux](https://github.com/tmux/tmux/wiki) or
[screen](https://www.gnu.org/software/screen/manual/screen.html#Overview) so
that you can run your workflow in a session in the background while doing
other things on the cluster or even logging out of the cluster.

## Run your workflow in an interactive job

For short workflows with only a few rules that need the same compute resources
in terms of CPU (cores), you can start an interactive job (in your `tmux` or
`screen` session) and run your Snakemake workflow as you would do that on your
local machine. Make sure to give your interactive job enough time to finish
running all rules of your Snakemake workflow.

## Use built-in SLURM support

For workflows with long run times and/or where each rule requires different
compute resources, Snakemake can be configured to automatically send each rule
as a job to the SLURM queue and to track the status of each job.

Since version 7.19.0 Snakemake comes with built-in support for execution on
compute clusters with the SLURM workload manager. To enable this you supply
the `--slurm` flag to your Snakemake command. In addition you need to
specify the id (_e.g. `snic-2023-01-001`) for your compute project. This can
be done directly on the command line with `--default-resources
slurm_account=snic2023-01-001`. You also need to specify the number of jobs
that Snakemake will queue at the same time with `-j`, _e.g._ `-j 100` to
allow up to 100 jobs to be put into the queue at any given time. So the
command would be (in addition to any other flags you may want to use):

```bash
snakemake --slurm --default-resources slurm_account=snic2023-01-001 -j 100
```

Snakemake will submit each job to the SLURM queue and inform you about both
the local jobid and the SLURM jobid by writing something similar to this to
your terminal:

```
Job 0 has been submitted with SLURM jobid 37099380 (log: .snakemake/slurm_logs/rule_name/37099380.log).
```

In this example the log output from the job will be in
`.snakemake/slurm_logs/rule_name/37099380.log`.

So how do you specify SLURM resources such as runtime, CPUs *etc.*? The best
way to do that is to use the `resources:` and `threads:` directives in the
rules of your workflow. This allows you to fine-tune jobs to run with
individual runtime and CPU usage. Take a look at the example rule below:

```python
rule testrule:
    output:
        "results/output.txt"
    resources:
        runtime = 60
    threads: 4
    shell:
        """
        uname -a > {output}
        """
```

This rule uses the standard resource `runtime` to set the maximum allowed
time for the rule to 60 minutes and sets the number of threads to 4. This
means that the rule will have a time limit of 60 minutes and will use 4 CPUs.

Of course you could set the runtime and threads using a configuration file
as we have shown in earlier sections of this tutorial, _e.g._ with a config
file that contains:

```yaml
testrule:
  threads: 4
  runtime: 60
```

```python
rule testrule:
    output:
        "results/output.txt"
    resources:
        runtime = config["testrule"]["runtime"]
    threads: config["testrule"]["threads"]
    shell:
        """
        uname -a > {output}
        """
```

Note that when using the `--slurm` flag `-j` only specifies the number of
jobs that can be sent to the queue at any given time, while the number of
CPUs used for each job is set via the `threads:` directive.

The `resources` directive can also be used to specify constraints, for
instance if jobs need to run on nodes with more memory you can use the
following on the Uppmax compute cluster:

```python
  resources:
    constraint = "mem256GB"
```

If you need to submit the job on another cluster, _e.g._ the 'snowy' cluster
on Uppmax you can do so with the `slurm_extra` keyword in the `resources`
directive:

```python
  resources:
    slurm_extra = "-M snowy"
```

You can read more details about running Snakemake on compute clusters in the
[Snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

## SLURM Profile

As an alternative to the built-in SLURM support you can also use a
configuration profile developed for SLURM, such as
[https://github.com/Snakemake-Profiles/slurm]() or the more light-weight
[https://github.com/jdblischak/smk-simple-slurm]().
