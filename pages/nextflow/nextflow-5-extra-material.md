The tutorial for Nextflow has intentionally been quite short, but there are many
more things you can do with it than what was covered here. In this section we
will briefly show some of these things if you are interested to know more
details about Nextflow, and we'll also provide some links to additional
resources below:

 * The Nextflow [documentation](https://www.nextflow.io/docs/latest/index.html)
 * [Learning Nextflow in 2020](https://www.nextflow.io/blog/2020/learning-nextflow-in-2020.html)
 * Common [Nextflow patterns](http://nextflow-io.github.io/patterns/index.html)
 * The [nf-core](https://nf-co.re/) pipeline collection
 * Nextflow training at [Seqera](https://seqera.io/training/)


## Running Nextflow on Uppmax

A lot of researchers in Sweden are using the Uppmax computer cluster in Uppsala,
which is easily handled by Nextflow. What you need to do is to add the following
*profile* to your `nextflow.config` file:

```
profiles {

    // Uppmax general profile
    uppmax {
        process {
            executor = 'slurm'
            clusterOptions = { '-A "account"' }
            memory = { 6.GB * task.attempt }
            cpus = { 1 * task.attempt }
            time = { 10.h * task.attempt }
            scratch = '$SNIC_TMP'
            errorStrategy  = 'retry'
            maxRetries = 1
        }
    }
}
```

This will add a profile to your workflow, which you can access by running the
workflow with `-profile uppmax`. You will have to edit the `"account"` part of
the code above to correspond to your project account, but the rest you can leave
as-is, unless you want to tinker with *e.g.* resource specifications. That's all
you need! Nextflow will take care of communications with SLURM (the system used
by Uppmax) and send off jobs to the cluster for you, and everything will look
exactly the same as if you were executing locally.

## Advanced channel creation

The input data shown in the MRSA example workflow is quite simple, but Nextflow
channels can do much more than that. A common scenario in high-throughput
sequencing is that you have pairs of reads for each sample. Nextflow has a
special, built-in way to create channels for this data type: the `fromFilePairs`
channel factory:

```groovy
Channel
    .fromFilePairs( "data/*_R{1,2}.fastq.gz" )
    .set( raw_reads )
```

This will create a channel containing all the reads in the `data/` directory on
the format `<sample>_R1.fastq.gz` or `<sample>_R2.fastq.gz` and pair them
together into a nested tuple looking like this:

```groovy
[sample, [data/sample_R1.fastq.gz, data/sample_R2.fastq.gz]]
```

The first index of the tuple (`[0]`) thus contains the value `sample`, while the
second index (`[1]`) contains another tuple with paths to both read files. This
nested tuple can easily be passed into processes for *e.g.* read alignment, and
it makes the entire procedure of going from read pairs (*i.e.* two separate
files, one sample) into a single alignment file (one file, one sample) very
simple.

We can also do quite advanced things when creating channels, such as this:

```groovy
Channel
    .fromPath( params.metadata )
    .splitCsv( sep: "\t", header: true )
    .map{ row -> tuple("${row.sample_id}", "${row.treatment}") }
    .filter{ it[1] != "DMSO" }
    .unique()
    .set { samples_and_treatments }
```

That's a bit of a handful! But what does it do? The first line specifies that we
want to read some data from a file specified by the `metadata` parameter, and
the second line actually reads that data using TAB as delimiter, including a
header. The `map` operator takes each entire row and subsets it to only two
columns: the `sample_id` and `treatment` columns. This subset is stored as a
tuple. The `filter` operator is then used to remove any tuples where the second
entry is equal to the string `"DMSO"` (*i.e.* untreated cells, in this example).
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
and send the resulting channels and their content to your downstream processes
in a simple way. While this toy example is probably too complicated for a lot of
projects, at least you know that there are many things you can do in regards to
input data complexity with Nextflow!

## Using Groovy in processes

You don't have to use bash or external scripts inside your processes all the
time unless you want to: Nextflow can use Groovy in the same way that Snakemake
uses Python. For example, look at this process:

```groovy
process index_fasta {
    tag "${fasta_name}"

    input:
    tuple val(fasta), path(fasta_file)

    output:
    path("${fasta_name}.idx")

    script:
    fasta_name = fasta.substring(0, fasta.lastIndexOf("."))
    """
    index --ref ${fasta_file},${fasta_name}
    """
}
```

Here we have some command `index` that, for whatever reason, requires both the
path to a FASTA file as well as the name of that file *without* the `.fasta`
extension. We can use Groovy in the `script` directive together with normal
bash: we can mix and match as we like. The first line of the `script` directive
gets the name of the FASTA file without the extension by removing anything after
the dot, while the second calls the `index` command like normal using Bash.
