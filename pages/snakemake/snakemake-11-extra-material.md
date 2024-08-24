If you want to read more about Snakemake in general you can find several
resources here:

- The Snakemake documentation is available on [ReadTheDocs](https://snakemake.readthedocs.io/en/stable/#).
- Here is another (quite in-depth) [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html#tutorial).
- If you have questions, check out [stack overflow](https://stackoverflow.com/questions/tagged/snakemake).

# Using containers in Snakemake

Snakemake also supports defining an Apptainer or Docker container for each rule
(you will have time to work on the [Containers
tutorial](containers-1-introduction) later during the course). Analogous to
using a rule-specific Conda environment, specify `container:
"docker://some-account/rule-specific-image"` in the rule definition. Instead of
a link to a container image, it is also possible to provide the path to a
`*.sif` file (= a _Singularity image file_). When executing Snakemake, add the
`--software-deployment-method apptainer` (or the shorthand `--sdm apptainer`)
flag to the command line. For the given rule,
an Apptainer container will then be created from the image or file that
is provided in the rule definition on the fly by Snakemake and the rule will be
run in this container.

You can find pre-made Apptainer or Docker images for many tools on
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
snakemake --software-deployment-method apptainer
```

Feel free to modify the MRSA workflow according to this example. As Apptainer
is a container software that was developed for HPC clusters, and for example the
Mac version is still a beta version, it might not work to run your updated
Snakemake workflow with Apptainer locally on your computer.
In the next section we explain how you can run Snakemake workflows on UPPMAX
where Apptainer is pre-installed.

# Running Snakemake workflows on HPC clusters

If you need to run a Snakemake workflow on a high-performance computing (HPC)
cluster you have a wide range of options at your disposal. Via the [plugin
catalog](https://snakemake.github.io/snakemake-plugin-catalog/) you can find
plugins that will add support for various HPC schedulers to Snakemake. 

Here we will focus on how to run Snakemake workflows on clusters with SLURM, a
workload manager commonly used on HPC clusters in Sweden such as
[Rackham](https://www.uu.se/centrum/uppmax/resurser/kluster/rackham),
[Tetralith](https://www.nsc.liu.se/systems/tetralith/) and
[Dardel](https://www.pdc.kth.se/hpc-services/computing-systems).

> **Tip**
> When running on remote clusters we highly recommend to use a session manager like
> [tmux](https://github.com/tmux/tmux/wiki) or
> [screen](https://www.gnu.org/software/screen/manual/screen.html#Overview) so
> that you can run your workflow in a session in the background while doing other
> things on the cluster or even logging out of the cluster.

## Option 1: Run the entire workflow as a single job

For short workflows with only a few rules that need the same compute resources
in terms of CPU (cores) and memory, you can submit the entire workflow as a job
directly to the SLURM scheduler, or start an interactive job (in your `tmux` or
`screen` session) and run your Snakemake workflow as you would do that on your
local machine. Make sure to give your job enough time to finish running all
rules of your Snakemake workflow.

If you choose this option, you don't need to install anything from the plugin
catalog. However, your workflow may not run as efficiently as it could if you
were to add SLURM support in Snakemake.

## Option 2: Use built-in SLURM support

For workflows with long run times and/or where each rule requires different
compute resources, Snakemake comes with built in functionality for interacting
with the SLURM workload manager and send each rule as a job to the SLURM queue
and to track the status of each job. 

In this case, you can start the workflow on the login node and let it run there
until all jobs have finished. Given that workflows often consist of many rules,
some of which may be highly resource demanding, this is the option we recommend
when running most Snakemake workflows on HPC clusters.

To add SLURM support to Snakemake you first need to install the [SLURM
plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
from the plugin catalog. This can be done with conda:

```bash
conda install -c conda-forge snakemake-executor-plugin-slurm
```

Once installed, adding the `--executor slurm` flag to your Snakemake command
line call will enable the plugin. You also need to specify how many jobs Snakemake
can submit to the SLURM queue at the same time with the `-j` flag. For example,
to allow up to 100 jobs to be put into the queue at any given time, you would
run Snakemake with the following command:

```bash
snakemake --executor slurm -j 100 <other flags>
```

### Specifying resources for SLURM

Depending on the cluster you are using, you will need to specify some resource
requirements for the rules in your workflow, such as the number of CPUs, memory,
runtime and account id. This can be done either:

1. directly on the command line with the `--default-resources` flag which sets
   default resource settings for all rules
2. in the rule definition of your workflow using the `resources:` directive, or
3. in a [configuration
   profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles),
   a folder with a `config.yaml` file that contains the resource settings.

You can also use a combination of these methods. For example, the SLURM account
id (_e.g. `naiss-2023-01-001`), which will most likely be the same for all
rules, can be set with `--default-resources`:

```bash
snakemake --executor slurm -j 100 --default-resources slurm_account=naiss-2023-01-001 
```

Rule-specific resources such as runtime, memory and number of CPUs can be set in
the rule definition, for example:

```python
rule testrule:
    output:
        "results/output.txt"
    resources:
        runtime = 60,
        mem_mb = 16000,
        cpus_per_task = 4
    shell:
        """
        uname -a > {output}
        """
```

This rule uses the standard resource `runtime` to set the maximum allowed time
(in minutes) for the rule, sets the memory requirement with `mem_mb` and the
number of requested CPUs with `cpus_per_task`. In this example the rule will have a
time limit of 60 minutes, will require 16G of RAM and 4 CPUs.

Some clusters also require you to specify the **partition** you want to run your
job on. The partition name will differ between clusters, for example the Rackham
cluster uses `core` and `node` partitions, while Dardel uses _e.g._ `shared` and
`main`. See the documentation for the cluster you are using for more
information.

The partition can be set with the `slurm_partition` resource, for example like so:

```python
rule testrule:
    output:
        "results/output.txt"
    resources:
        runtime = 60,
        mem_mb = 16000,
        cpus_per_task = 4,
        slurm_partition: "shared"
    shell:
        """
        uname -a > {output}
        """
```

To make it easy to adapt your workflow to different compute clusters it is
recommended to define resource settings in a **configuration profile**. A
configuration profile is a folder with a `config.yaml` file that contains values
for Snakemake command line arguments, allowing you to modify the behavior of
Snakemake without changing the workflow code. For example, you could create a
`dardel` folder (_e.g._ in the root of your workflow) with a `config.yaml` file
that contains the following:

```yaml
executor: "slurm"
jobs: 100
default-resources:
  slurm_account: "naiss-2023-01-001"
  slurm_partition: "shared"
  mem_mb: 16000
  cpus_per_task: 4
  runtime: 60
```

This yaml-formatted file contains Snakemake command line arguments that will be used
when running the workflow. You can then run Snakemake with the `--profile` flag pointing
to the folder containing the `config.yaml` file:

```bash
snakemake --profile dardel
```

This greatly simplifies running the workflow on different clusters, and makes
the command line call much more succinct.

To set rule-specific resources in the configuration profile, you can add a
`set_resources:` section to the `config.yaml` file:

```yaml
executor: "slurm"
jobs: 100
default-resources:
  slurm_account: "naiss-2023-01-001"
  slurm_partition: "shared"
  mem_mb: 16000
  cpus_per_task: 4
  runtime: 60
set_resources:
  index_genome:
    runtime: 240
    mem_mb: 32000
    cpus_per_task: 8
  align_to_genome:
    runtime: 120
    mem_mb: 24000
    cpus_per_task: 6
```

In this example, the `index_genome` rule will have a runtime of 240 minutes,
will require 32G of RAM and 8 CPUs, while the `align_to_genome` rule will have a
runtime of 120 minutes, will require 24G of RAM and 6 CPUs. Both rules will use
the `slurm_account` and `slurm_partition` settings from the `default_resources`
section, unless overridden in the rule-specific settings.

You can still define resources in the rule definition, but the values in the
configuration profile will take precedence.

Now, when you run your Snakemake workflow with:

```bash
snakemake --profile dardel
```

Snakemake will submit each job to the SLURM queue and inform you about both
the local jobid and the SLURM jobid by writing something similar to this to
your terminal:

```
Job 0 has been submitted with SLURM jobid 37099380 (log: .snakemake/slurm_logs/rule_name/37099380.log).
```

In this example the log output from the job will be in
`.snakemake/slurm_logs/rule_name/37099380.log`.

You can read more details about running Snakemake on compute clusters in the
[Snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).