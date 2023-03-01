It's generally a good idea to separate project-specific parameters from the
actual implementation of the workflow. If we want to move all project-specific
information to `config.yml`, and let the Snakefile be a more general RNA-seq
analysis workflow, we need the config file to:

* Specify which samples to run.
* Specify which genome to align to and where to download its sequence and
  annotation files.
* (Any other parameters we might need to make it into a general workflow,
  *e.g.* to support both paired-end and single-read sequencing)

!!! Note
    Putting all configuration in `config.yml` will break the
    `generate_rulegraph` rule. You can fix it either by replacing
    `--config max_reads=0` with `--configfile=config.yml` in the shell
    command of that rule in the Snakefile, or by adding
    `configfile: "config.yml"` to the top of the Snakefile (as mentioned
    in a previous tip).

The first point is straightforward; rather than using `SAMPLES = ["..."]` in
the Snakefile we define it as a parameter in `config.yml`. You can either add
it as a list similar to the way it was expressed before by adding
 `SAMPLES: ["..."]` to `config.yml`, or you can use this yaml notation:

```yaml
sample_ids:
- SRR935090
- SRR935091
- SRR935092
```

You also have to change the workflow to reference `config["sample_ids"]` (if
using the latter example) instead of `SAMPLES`, as in:

```bash
expand("intermediate/{sample_id}_fastqc.zip",
            sample_id = config["sample_ids"])
```

Do a dry-run afterwards to make sure that everything works as expected.

The second point is trickier. Writing workflows in Snakemake is quite
straightforward when the logic of the workflow is reflected in the file names,
*i.e.* `my_sample.trimmed.deduplicated.sorted.fastq`, but that isn't always the
case. In our case we have the FTP paths to the genome sequence and annotation
where the naming doesn't quite fit with the rest of the workflow. The easiest
solution is probably to make three parameters to hold these values, say
`genome_id`, `genome_fasta_path` and `genome_gff_path`, but we will go for
a somewhat more complex but very useful alternative. We want to construct
a dictionary where something that will be a wildcard in the workflow is the key
and the troublesome name is the value. An example might make this clearer (this
is also in `config.yml` in the finished version of the workflow under
`tutorials/git/`). This is a nested dictionary where "genomes" is a key with
another dictionary as value, which in turn has genome ids as keys and so on. The
idea is that we have a wildcard in the workflow that takes the id of a genome as
value (either "NCTC8325" or "ST398" in this case). The fasta and gff3 paths can
then be retrieved based on the value of the wildcard.

```yaml
genomes:
  NCTC8325:
    fasta: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
    gff3: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz
  ST398:
    fasta: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection//staphylococcus_aureus_subsp_aureus_st398/dna/Staphylococcus_aureus_subsp_aureus_st398.ASM958v1.dna.toplevel.fa.gz
    gff3: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_st398//Staphylococcus_aureus_subsp_aureus_st398.ASM958v1.37.gff3.gz
```

Go ahead and add the section above to `config.yml`.

Let's now look at how to do the mapping from genome id to fasta path in the
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
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O {output} -o {log}
        """
```

We don't want the hardcoded genome id `NCTC8325`, so replace that with a
wildcard, say `{genome_id}` (remember to add the wildcard to the `log:`
directive as well).

Also change in `index_genome` to use a wildcard rather than a hardcoded genome
id. Here you will run into a complication if you have followed the previous
instructions and use the `expand()` expression. We want the list to expand to
`["intermediate/{genome_id}.1.bt2", "intermediate/{genome_id}.2.bt2", ...]`,
*i.e.* only expanding the wildcard referring to the bowtie2 index. To keep the
`genome_id` wildcard from being expanded we have to "mask" it with double curly
brackets: `{{genome_id}}`. In addition, we need to replace the hardcoded
`intermediate/NCTC8325` in the shell directive of the rule with the genome id
wildcard. Inside the shell directive the wildcard object is accessed with this
syntax: `{wildcards.genome_id}`, so the bowtie2-build command should be:

```bash
bowtie2-build tempfile intermediate/{wildcards.genome_id} > {log}
```

We now need to supply the remote paths to the fasta and gff files for a given
genome id. Because we've added this information to the config file we just need
to pass it to the rule in some way.

Take a look at the code and `get_genome_fasta` rule below. Here we have defined
a function called `get_fasta_path` which takes the `wildcards` object as its
only argument. This object allows access to the wildcards values via attributes
(here `wildcards.genome_id`). The function will then look in the nested `config`
dictionary and return the value of the fasta path for the key
`wildcards.genome_id`. In the rule this path is stored in the `fasta_path` param
value and is made available to `wget` in the shell directive.

```python
def get_fasta_path(wildcards):
    return config["genomes"][wildcards.genome_id]["fasta"]

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/{genome_id}.fa.gz"
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    params:
        fasta_path = get_fasta_path
    shell:
        """
        wget {params.fasta_path} -O {output} -o {log}
        """
```

Note that this will only work if the `{genome_id}` wildcard can be resolved to
something defined in the config (currently `NCTC8325` or `ST398`). If you try to
generate a fasta file for a genome id not defined in the config Snakemake will
complain, even at the dry-run stage.

Now change the `get_genome_gff3` rule in a similar manner.

The rules `get_genome_fasta`, `get_genome_gff3` and `index_genome` can now
download and index *any genome* as long as we provide valid links in the config
file.

However, we need to define somewhere which genome id we actually want to use
when running the workflow. This needs to be done both in `align_to_genome` and
`generate_count_table`. Do this by introducing a parameter in `config.yml`
called `"genome_id"` (you can set it to either `NCTC8325` or `ST398`).

Now we can resolve the `genome_id` wildcard from the config. See below for an
example for `align_to_genome`. Here the `substr` wildcard gets expanded from a
list while `genome_id` gets expanded from the config file.

```python
input:
    index = expand("intermediate/{genome_id}.{substr}.bt2",
           genome_id = config["genome_id"],
           substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
```

Also change the hardcoded genome id in the `generate_count_table` input in a
similar manner.

In general, we want the rules as far downstream as possible in the workflow to
be the ones that determine what the wildcards should resolve to. In our case
this is `align_to_genome` and `generate_count_table`. You can think of it like
the rule that really "needs" the file asks for it, and then it's up to Snakemake
to determine how it can use all the available rules to generate it. Here the
`align_to_genome` rule says "I need this genome index to align my sample to" and
then it's up to Snakemake to determine how to download and build the index.

One last thing is to change the hardcoded `NCTC8325` in the `shell:` directive
of `align_to_genome`. Bowtie2 expects the index name supplied with the `-x` flag
to be without the ".*.bt2" suffix so we can't use `-x {input.index}`. Instead
we'll insert the genome_id directly from the config like this:

```bash
shell:
    """  
    bowtie2 -x intermediate/{config[genome_id]} -U {input[0]} > {output} 2>{log}
    """
```

!!! Success "Quick recap"
    In this section we've learned how to generalize a Snakemake workflow with a
    number of excellent features:

    - A general RNA-seq pipeline which can easily be reused between projects,
    thanks to clear separation between code and settings.
    - Great traceability due to logs and summary tables.
    - Clearly defined the environment for the workflow using Conda.
    - The workflow is neat and free from temporary files due to using `temp()` and
    `shadow`.
    - A logical directory structure which makes it easy to separate raw data,
    intermediate files, and results.
    - A project set up in a way that makes it very easy to distribute and
    reproduce either via Git, Snakemake's `--archive` option or a Docker image.
