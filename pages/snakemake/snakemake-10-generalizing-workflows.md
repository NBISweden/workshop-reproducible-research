It's a good idea to separate project-specific parameters from the
actual implementation of the workflow. This allows anyone using the
workflow to modify its behaviour without changing the underlying code,
making the workflow more general.

In order to generalize our RNA-seq analysis workflow we should move all
project-specific information to `config.yml`. This means that we want the
config file to:

* Specify which samples to run.
* Specify which genome to align to and where to download its sequence and
  annotation files.
* (Contain any other parameters we might need to make it into a general
  workflow, *e.g.* to support both paired-end and single-read sequencing)

> **Note** <br>
> Putting all configuration in `config.yml` will break the
> `generate_rulegraph` rule. You can fix it either by replacing
> `--config max_reads=0` with `--configfile=config.yml` in the shell
> command of that rule in the Snakefile, or by adding
> `configfile: "config.yml"` to the top of the Snakefile (as mentioned
> in a previous tip).

The first point is straightforward; rather than using `SAMPLES = ["..."]` in
the Snakefile we define it as a parameter in `config.yml`. You can either add
it as a list similar to the way it was expressed before by adding:

```yaml
SAMPLES: ["SRR935090", "SRR935091", "SRR935092"]
```

To `config.yml`, or you can use this YAML notation (whether you
choose `SAMPLES` or `sample_ids` as the name of the entry doesn't matter,
you will just have to reference the same name in the config dictionary
inside the workflow):

```yaml
sample_ids:
  - SRR935090
  - SRR935091
  - SRR935092
```

Change the workflow to reference `config["sample_ids"]` (if using the latter
example) instead of `SAMPLES`, as in:

```bash
expand("results/fastqc/{sample_id}_fastqc.zip",
            sample_id = config["sample_ids"])
```

Remove the line with `SAMPLES = ["SRR935090", "SRR935091", "SRR935092"]`
that we added to the top of `snakefile_mrsa.smk` in
[Snakemake 8: Targets](snakemake-8-targets).

Do a dry-run afterwards to make sure that everything works as expected.

You may remember from the [snakemake-5-parameters](snakemake-5-parameters)
part of this tutorial that we're using a function to return the URL of the
FASTQ files to download for each sample:

```python
def get_sample_url(wildcards):
    samples = {
        "SRR935090": "https://figshare.scilifelab.se/ndownloader/files/39539767",
        "SRR935091": "https://figshare.scilifelab.se/ndownloader/files/39539770",
        "SRR935092": "https://figshare.scilifelab.se/ndownloader/files/39539773"
    }
    return samples[wildcards.sample_id]
```

Here the URLs of each sample_id is hard-coded in the `samples` dictionary
inside the function. To generalize this function we can move the definition
to the config file, placing it for example under an entry that we call
`sample_urls` like this:

```yaml
sample_urls:
  SRR935090: "https://figshare.scilifelab.se/ndownloader/files/39539767"
  SRR935091: "https://figshare.scilifelab.se/ndownloader/files/39539770"
  SRR935092: "https://figshare.scilifelab.se/ndownloader/files/39539773"
```

This is what's called 'nested' key/value pairs, meaning that each sample_id
-> URL pair becomes nested under the config key `sample_urls`. So in order
to access the URL of _e.g._ `SRR935090` we would use
`config["sample_urls"]["SRR935090"]`. This means that you will have to
update the `get_sample_url` function to:

```python
def get_sample_url(wildcards):
    return config["sample_urls"][wildcards.sample_id]
```

Now the function uses the global `config` dictionary to return URLs for each
sample_id. Again, do a dry-run to see that the new implementation works.

> *Tip!* <br>
> If you were to scale up this workflow with more samples it could become
> impractical to have to define the URLs by hand in the config file. A
> tip then is to have a separate file where samples are listed in one column
> and the URLs (or file paths) in another column. With a few lines of python
> code you could then read that list at the start of the workflow and add
> each sample to the config dictionary.

Now let's take a look at the genome reference used in the workflow. In the
`get_genome_fasta` and `get_genome_gff3` rules we have hard-coded FTP paths
to the FASTA GFF annotation file for the genome `NCTC8325`. We can
generalize this in a similar fashion to what we did with the
`get_SRA_by_accession` rule. Let's add a nested entry called `genomes` to
the config file that will hold the genome id and FTP paths to the FASTA and
GFF file:

```yaml
genomes:
  NCTC8325:
    fasta: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
    gff3: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz
  ST398:
    fasta: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection//staphylococcus_aureus_subsp_aureus_st398/dna/Staphylococcus_aureus_subsp_aureus_st398.ASM958v1.dna.toplevel.fa.gz
    gff3: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_st398//Staphylococcus_aureus_subsp_aureus_st398.ASM958v1.37.gff3.gz
```

As you can see this is very similar to what with did with `sample_urls`,
just that we have one more nested level. Now to access the FTP path to the
FASTA file for genome id `NCTC8325` we can use
`config["genomes"]["NCTC8325"]["fasta"]`.

Let's now look at how to do the mapping from genome id to FASTA path in the
rule `get_genome_fasta`. This is how the rule currently looks (if you have
added the log section as previously described).

```python
rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/NCTC8325.fa.gz"
    log:
        "results/logs/get_genome_fasta/NCTC8325.log"
    shell:
        """
        wget -o {log} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O {output}
        """
```

We don't want the hard-coded genome id `NCTC8325`, so replace that with a
wildcard, say `{genome_id}` (remember to add the wildcard to the `log:`
directive as well). We now need to supply the remote paths to the FASTA file
for a given genome id. Because we've added this information to the
config file we just need to pass it to the rule in some way, and just like
in the `get_SRA_by_accession` rule we'll use a function to do the job:


```python
def get_fasta_path(wildcards):
    return config["genomes"][wildcards.genome_id]["fasta"]

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/ref/{genome_id}.fa.gz"
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    params:
        fasta_path = get_fasta_path
    shell:
        """
        wget -o {log} {params.fasta_path} -O {output}
        """
```

Now change the `get_genome_gff3` rule in a similar manner. Click to see the
solution below if you're having trouble.

<details>
<summary> "Click to see solution" </summary>

```python
def get_gff_path(wildcards):
    return config["genomes"][wildcards.genome_id]["gff3"]

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/ref/{genome_id}.gff3.gz"
    log:
        "results/logs/get_genome_gff3/{genome_id}.log"
    params:
        gff3_path = get_gff_path
    shell:
        """
        wget -o {log} {params.gff3_path} -O {output}
        """
```
</details>

Also change in `index_genome` to use a wildcard rather than a hard-coded genome
id. Here you will run into a complication if you have followed the previous
instructions and use the `expand()` expression. We want the list to expand to
`["results/bowtie2/{genome_id}.1.bt2", "results/bowtie2/{genome_id}.2.bt2", ...]`,
*i.e.* only expanding the wildcard referring to the Bowtie2 index. To keep the
`genome_id` wildcard from being expanded we have to "mask" it with double curly
brackets: `{{genome_id}}`. In addition, we need to replace the hard-coded
`results/bowtie2/NCTC8325` in the shell directive of the rule with the genome id
wildcard. Inside the shell directive the wildcard object is accessed with this
syntax: `{wildcards.genome_id}`, so the Bowtie2-build command should be:

```bash
bowtie2-build tempfile results/bowtie2/{wildcards.genome_id} > {log}
```

Note that this will only work if the `{genome_id}` wildcard can be resolved to
something defined in the config (currently `NCTC8325` or `ST398`). If you try to
generate a FASTA file for a genome id not defined in the config Snakemake will
complain, even at the dry-run stage.

The rules `get_genome_fasta`, `get_genome_gff3` and `index_genome` can now
download and index *any genome* as long as we provide valid links in the config
file.

However, we need to define somewhere which genome id we actually want to use
when running the workflow. This needs to be done both in `align_to_genome` and
`generate_count_table`. Do this by introducing a parameter in `config.yml`
called `"genome_id"` (you can set it to either `NCTC8325` or `ST398`), _e.g._:

```yaml
genome_id: "NCTC8325"
```

Now we can resolve the `genome_id` wildcard from the config. See below for an
example for `align_to_genome`. Here the `substr` wildcard gets expanded from a
list while `genome_id` gets expanded from the config file.

```python
input:
    index = expand("results/bowtie2/{genome_id}.{substr}.bt2",
           genome_id = config["genome_id"],
           substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
```

Also change the hard-coded genome id in the `generate_count_table` input in a
similar manner.

In general, we want the rules as far downstream as possible in the workflow to
be the ones that determine what the wildcards should resolve to. In our case
this is `align_to_genome` and `generate_count_table`. You can think of it like
the rule that really "needs" the file asks for it, and then it's up to Snakemake
to determine how it can use all the available rules to generate it. Here the
`align_to_genome` rule says "I need this genome index to align my sample to" and
then it's up to Snakemake to determine how to download and build the index.

One last thing is to change the hard-coded `NCTC8325` in the `shell:` directive
of `align_to_genome`. Bowtie2 expects the index name supplied with the `-x` flag
to be without the ".*.bt2" suffix so we can't use `-x {input.index}`. Instead
we'll insert the genome_id directly from the config like this:

```bash
shell:
    """
    bowtie2 -x results/bowtie2/{config[genome_id]} -U {input[0]} > {output} 2>{log}
    """
```

> **Summary** <br>
> Well done! You now have a complete Snakemake workflow with a number of
> excellent features:
>
> - A general RNA-seq pipeline which can easily be reused between projects,
>   thanks to clear separation between code and settings.
> - Great traceability due to logs and summary tables.
> - Clearly defined the environment for the workflow using Conda.
> - The workflow is neat and free from temporary files due to using `temp()` and
>   `shadow`.
> - A logical directory structure which makes it easy to separate data and
>   results of different software packages.
> - A project set up in a way that makes it very easy to distribute and
>   reproduce either via Git, Snakemake's `--archive` option or a Docker image.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to generalize a Snakemake workflow.
