# Extra material

## Running jobs in Singularity or Docker containers

Snakemake also supports defining a Singularity or Docker container for each rule
(you will have time to work on the [Docker tutorial](docker.md) and the
[Singularity tutorial](singularity.md) later during the course).
Analogous to using a rule-specific Conda environment, specify
`container: "docker://some-account/rule-specific-image"` in the rule definition.
Instead of a link to a container image, it is also possible to provide the path
to a `*.sif` file (= a Singularity file).
When executing Snakemake, add the `--use-singularity` flag to the command line.
For the given rule, a Singularity container will then be created from the image
or Singularity file that is provided in the rule definition on the fly by Snakemake
and the rule will be run in this container.

You can find pre-made Singularity or Docker images for many tools on [https://biocontainers.pro/](https://biocontainers.pro/)
(bioinformatics-specific) or on [https://hub.docker.com/](https://hub.docker.com/).

Here is an example for a rule and its execution:

```python
rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    input:
        "data/raw_internal/{sra_id}.fastq.gz",
        "intermediate/NCTC8325.1.bt2",
        "intermediate/NCTC8325.2.bt2",
        "intermediate/NCTC8325.3.bt2",
        "intermediate/NCTC8325.4.bt2",
        "intermediate/NCTC8325.rev.1.bt2",
        "intermediate/NCTC8325.rev.2.bt2"
    output:
        "intermediate/{sra_id,\w+}.bam"
    container: "docker://quay.io/biocontainers/bowtie2:2.3.4.1--py35h2d50403_1"
    shell:
        """
        bowtie2 -x intermediate/NCTC8325 -U {input[0]} > {output}
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

## Running Snakemake workflows on UPPMAX

There are several options to execute Snakemake workflows on UPPMAX (a HPC
cluster with the SLURM workload manager). In any case, we highly recommend to use
a session manager like `tmux` or `screen` so that you can run your workflow in a
session in the background while doing other things on the cluster or even logging
out of the cluster.

### Run your workflow in an interactive job

For short workflows with only a few rules that need the same compute resources
in terms of CPU (cores), you can start an interactive job (in your `tmux` or
`screen` session) and run your Snakemake workflow as you would do that on your
local machine. Make sure to give your interactive job enough time to finish
running all rules of your Snakemake workflow.

### Cluster configuration

For workflows with long run times and/or where each rule requires different
compute resources, Snakemake can be configured to automatically send each rule
as a job to the SLURM queue and to track the status of each job.

The relevant parameters for such a cluster configuration are `--cluster` and
`--cluster-config`, in combination with a `cluster.yaml` file that specifies
default and rule-specific compute resources and your compute account details.

Here is an example for a `cluster.yaml` file:

```yaml
# cluster.yaml - cluster configuration file
__default__:
  account: # fill in your project compute account ID
  partition: core
  time: 01:00:00
  ntasks: 1
  cpus-per-task: 1
### rule-specific resources
trimming:
  time: 01-00:00:00
mapping:
  time: 01-00:00:00
  cpus-per-task: 16
```

Start your Snakemake workflow in a `tmux` or `screen` session with the
following command:

```bash
snakemake
    -j 10 \
    --cluster-config cluster.yaml \
    --cluster "sbatch \
               -A {cluster.account} \
               -p {cluster.partition} \
               -t {cluster.time} \
               --ntasks {cluster.ntasks} \
               --cpus-per-task {cluster.cpus-per-task}"
```

The additional parameter `-j` specifies the number of jobs that Snakemake is
allowed to send to SLURM at the same time.

### SLURM Profile

In future Snakemake versions, the cluster configuration will be replaced
by so-called Profiles. The SLURM Profile needs to be set up with the software
[cookiecutter](https://cookiecutter.readthedocs.io/en/1.7.2/) (available via Conda).
During the [SLURM Profile](https://github.com/Snakemake-Profiles/slurm) setup,
you will be asked for several values for your Profile, e.g. for a Profile name
or your compute project account ID.

Rule-specific resources can be defined in each rule via the `resources: `
directive, for example:

```python
rule align_to_genome:
    input:
        "{genome_id}.bt2",
        "{sample}.fastq.gz"
    output:
        "{sample}.bam"
    resources:
        runtime = 360
    threads: 10
    shell:
        """
        aligner -t {threads} -i {input[1]} -x {input[0]} > {output}
        """
```

Any rule for which runtime is specified in the `resources` directive will be
submitted as one job to the SLURM queue with runtime as the allocated time.
Similarly, the number specified in the `threads` directive will be used as the
number of allocated cores.

You can even specify if rules should be resubmitted to the SLURM queue,
asking for more resources on subsequent attempts. Do this by modifying the
`resources` directive with _e.g._:

```python
    resources:
        runtime = lambda wildcards, attempt: attempt*360
```

In this example, the rule will ask for 360 minutes on the first attempt,
for 2*360 minutes on the second attempt, etc.

Finally, instead of specifying compute resources in the `resources` and
`threads` directives, it is possible to set up the SLURM Profile by providing
a `cluster.yaml` file. When you create the profile with cookiecutter and you
are prompted for `cluster_config []:` enter the path to a `cluster.yaml` file,
_e.g._: `cluster_config []: config/cluster.yaml`. Now you can control resource
allocations for rules either using the `resources:` directive in each rule,
_or_ by adding that information to the `config/cluster.yaml` file. If resources
are found in both locations, the allocations in the `cluster.yaml` file will
take precedence.

With this setup you can start the workflow with your SLURM Profile as follows
from within a `tmux` or `screen` session:

```bash
snakemake -j 10 --profile your_profile_name
```
