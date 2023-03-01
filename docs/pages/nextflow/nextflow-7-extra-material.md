There are many more things you can do with Nextflow than covered
here. If you are interested to learn more details about Nextflow, we will
briefly show some of its advanced features in this section. But first, here are
some links to additional resources on Nextflow:

 * [Nextflow patterns](http://nextflow-io.github.io/patterns/index.html) that
   can help with common operations and concepts
 * The Nextflow [documentation](https://www.nextflow.io/docs/latest/index.html)
 * [Learning Nextflow in 2020](https://www.nextflow.io/blog/2020/learning-nextflow-in-2020.html)
 * Nextflow training at [Seqera](https://seqera.io/training/)
 * A work-in-progress [Nextflow Carpentry course](https://carpentries-incubator.github.io/workflows-nextflow/index.html)
 * Community help from [Nextflow's Slack channel](https://join.slack.com/t/nextflow/shared_invite/zt-11iwlxtw5-R6SNBpVksOJAx5sPOXNrZg)

# Using containers in Nextflow

Nextflow has built-in support for using both Docker and Singularity containers (and others too),
either with a single container for the workflow as a whole or separate
containers for each individual process. The simplest way to do it is to have a
single container for your entire workflow, in which case you simply run the
workflow and specify the image you want to use, like so:

```bash
# Run with docker
nextflow run main.nf -with-docker [image]

# Run with Singularity
nextflow run main.nf -with-singularity [image].sif
```

If you don't want to supply this at every execution, you can also add it
directly to your configuration file:

```bash
# Docker configuration
process.container = 'image'
docker.enabled = true

# Singularity configuration
process.container = 'path/to/image.sif'
singularity.enabled = true
```

If you instead would like to have each process use a different container you can
use the `container` directive in your processes:

```groovy
process PROCESS_01 {
    (...)
    container: 'image_01'
    (...)
}

process PROCESS_02 {
    (...)
    container: 'image_02'
    (...)
}
```

Regardless of which solution you go for, Nextflow will execute all the processes
inside the specified container. In practice, this means that Nextflow will
automatically wrap your processes and run them by executing the Docker or
Singularity command with the image you have provided.

# Using Conda in Nextflow

While you can execute Nextflow inside Conda environments just like you would any
other type of software, you can also use Conda with Nextflow in the same way as
for Docker and Singularity above. You can either supply an `environment.yml`
file, the path to an existing environment or the packages and their versions
directly in the `conda` directive, like so:

```groovy
process PROCESS_01 {
    (...)
    conda: 'mrsa-environment.yml'
    (...)
}
process PROCESS_02 {
    (...)
    conda: 'path/to/mrsa-env'
    (...)
}
process PROCESS_03 {
    (...)
    conda: 'bioconda::bwa=0.7.17 bioconda::samtools=1.13'
    (...)
}
```

You can use either of the methods described above with your configuration file
as well, here exemplified using an `environment.yml` file:

```groovy
process.conda = 'mrsa-environment.yml'
```

# Running Nextflow on Uppmax

A lot of researchers in Sweden are using the Uppmax computer cluster in Uppsala,
which is easily handled by Nextflow. What you need to do is to add the following
*profile* to your `nextflow.config` file:

```
profiles {

    // Uppmax general profile
    uppmax {
        params{
            account        = null
        }
        process {
            executor       = 'slurm'
            clusterOptions = "-A '${params.account}'"
            memory         = { 6.GB * task.attempt }
            cpus           = { 1 * task.attempt }
            time           = { 10.h * task.attempt }
            scratch        = '$SNIC_TMP'
            errorStrategy  = 'retry'
            maxRetries     = 1
        }
    }
}
```

This will add a profile to your workflow, which you can access by running the
workflow with `-profile uppmax`. You will also have to supply an extra parameter 
`account` which corresponds to your SNIC project account, but the rest you can
leave as-is, unless you want to tinker with *e.g.* compute resource
specifications. That's all you need! Nextflow will take care of communications
with SLURM (the system used by Uppmax, specified by the `executor` line) and
will send off jobs to the cluster for you, and everything will look exactly the
same way as if you were executing the pipeline locally.

The `memory`, `cpus` and `time` lines define the various resources Nextflow will
use as well as how much to automatically increase them by if re-trying failed
tasks; this, in turn, is specified by the `errorStrategy` and `maxRetries`
variables. The `scratch` variable defines where each node's local storage is
situated, which gives Nextflow the most optimal access to the Uppmax file system
for temporary files.

# Advanced channel creation

The input data shown in the MRSA example workflow is not that complex, but
Nextflow channels can do much more than that. A common scenario in
high-throughput sequencing is that you have pairs of reads for each sample.
Nextflow has a special, built-in way to create channels for this data type: the
`fromFilePairs` channel factory:

```groovy
Channel
    .fromFilePairs ( "data/*_R{1,2}.fastq.gz" )
    .set           { ch_raw_reads }
```

This will create a channel containing all the reads in the `data/` directory in
the format `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` and will pair them
together into a nested tuple looking like this:

```groovy
[sample, [data/sample_R1.fastq.gz, data/sample_R2.fastq.gz]]
```

The first element of the tuple (`[0]`) thus contains the value `sample`, while
the second element (`[1]`) contains another tuple with paths to both read files.
This nested tuple can be passed into processes for *e.g.* read alignment, and
it makes the entire procedure of going from read pairs (*i.e.* two separate
files, one sample) into a single alignment file (one file, one sample) very
simple. For more methods of reading in data see the Nextflow documentation on
[Channel Factories](https://www.nextflow.io/docs/latest/channel.html#channel-factory).

We can also do quite advanced things to manipuate data in channels, such as this:

```groovy
Channel
    .fromPath ( params.metadata )
    .splitCsv ( sep: "\t", header: true )
    .map      { row -> tuple("${row.sample_id}", "${row.treatment}") }
    .filter   { id, treatment -> treatment != "DMSO" }
    .unique   (  )
    .set      { samples_and_treatments }
```

That's a bit of a handful! But what does it do? The first line specifies that we
want to read some data from a file specified by the `metadata` parameter, and
the second line actually reads that data using tab as delimiter, including a
header. The `map` operator takes each entire row and subsets it to only two
columns: the `sample_id` and `treatment` columns. This subset is stored as a
tuple. The `filter` operator is then used to remove any tuples where the second
entry, `treatment`, is not equal to the string `"DMSO"` (*i.e.* untreated cells, in this example).
We then only take the unique tuples and set the results as the new channel
`samples_and_treatments`. Let's say that this is the metadata we're reading:

```no-highlight
sample_id     dose    group     treatment
sample_1      0.1     control   DMSO
sample_1      1.0     control   DMSO
sample_1      2.0     control   DMSO
sample_2      0.1     case      vorinostat
sample_2      1.0     case      vorinostat
sample_2      2.0     case      vorinostat 
sample_3      0.1     case      fulvestrant
sample_3      1.0     case      fulvestrant
sample_3      2.0     case      fulvestrant
```

Given the channel creation strategy above, we would get the following result:

```no-highlight
[sample_2, vorinostat]
[sample_3, fulvestrant]
```

In this way, you can perform complex operations on input files or input metadata
and send the resulting content to your downstream processes
in a simple way. Composing data manipuations in Nextflow like this can be half 
the fun of writing the workflow. Check out Nextflow's documentation on 
[Channel operators](https://www.nextflow.io/docs/latest/operator.html) to
see the full list of channel operations at your disposal. 

# Using Groovy in processes

You don't have to use bash or external scripts inside your processes all the
time unless you want to: Nextflow is based on Groovy, which allows you to use
both Groovy and Bash in the same process. For example, have a look at this:

```groovy
process index_fasta {
    tag "${fasta_name}"

    input:
    tuple val(fasta), path(fasta_file)

    output:
    path("${fasta_name}.idx"), emit: fasta

    script:
    fasta_name = fasta.substring(0, fasta.lastIndexOf("."))
    """
    index --ref ${fasta_file},${fasta_name}
    """
}
```

Here we have some command `index` that, for whatever reason, requires both the
path to a FASTA file and the name of that file *without* the `.fasta`
extension. We can use Groovy in the `script` directive together with normal
Bash, mixing and matching as we like. The first line of the `script` directive
gets the name of the FASTA file without the extension by removing anything
after the dot, while the second calls the `index` command like normal using
bash.

# The nf-core pipeline collection

You may have heard of the [nf-core](https://nf-co.re/) pipeline collection
previously, which is a large, collaborative bioinformatics community dedicated
to building, developing and maintaining Nextflow workflows. In fact, if you have
sequenced data at *e.g.* the National Genomics Infrastructure ([NGI](https://ngisweden.scilifelab.se/)),
you can be sure that the data processing has been run using one of the nf-core
pipelines! While the community only started in 2018 (with a [Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x)
paper in 2020), it already has over 30 production-ready pipelines with
everything from genomics, transcriptomics, proteomics and metagenomics - and
more being developed all the time.

The nf-core pipelines all work in the same way, in that they have the same exact
base for inputs, parameters and arguments, making them all highly similar to
run. Since you've already learnt the basics of Nextflow in this course, you
should now be able to also run the nf-core pipelines! It might be that you have
a data type that you can analyse using one of the pipelines in nf-core, meaning
you don't need to do anything other than find out what parameters you should run
it with.

Each pipeline comes with extensive documentation, test datasets that you can
use to practice on, can be run on both HPCs like Uppmax, cloud services like
AWS or locally on your own computer. All pipelines support both Conda and
Docker/Singularity, and you can additionally run specific versions of the
pipelines, allowing for full reproducibility of your analyses. If you want to
check nf-core out, simply head over to their [list of pipelines](https://nf-co.re/pipelines)
and see what's available! Who knows, you might even write your own nf-core
pipeline in the future?
