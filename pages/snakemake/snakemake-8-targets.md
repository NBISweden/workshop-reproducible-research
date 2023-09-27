We've mentioned that Snakemake rules take either strings or a list of strings as
input, and that we can use any Python expression in Snakemake workflows. 
Here we'll show how these features help us condense the code of rules.

Consider the rule `align_to_genome` below.

```python
rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        "intermediate/{sample_id}.bam"
    input:
        "data/raw_internal/{sample_id}.fastq.gz",
        "intermediate/NCTC8325.1.bt2",
        "intermediate/NCTC8325.2.bt2",
        "intermediate/NCTC8325.3.bt2",
        "intermediate/NCTC8325.4.bt2",
        "intermediate/NCTC8325.rev.1.bt2",
        "intermediate/NCTC8325.rev.2.bt2"
    shell:
        """
        bowtie2 -x intermediate/NCTC8325 -U {input[0]} > {output}
        """
```

Here we have seven inputs; the FASTQ file with the reads and six files with
similar file names from the Bowtie2 genome indexing. Instead of writing all 
the filenames we can tidy this up by using a Python expression to generate a 
list of these files instead. If you're familiar with Python you could do 
this with 
[list comprehensions](https://docs.python.org/3/tutorial/datastructures.html#list-comprehensions)
like this:

```python
input:
    "data/raw_internal/{sample_id}.fastq.gz",
    [f"intermediate/NCTC8325.{substr}.bt2" for
        substr in ["1", "2", "3", "4", "rev.1", "rev.2"]]
```

This will take the elements of the list of substrings one by one, and insert
that element in the place of `{substr}`. Since this type of aggregating
rules are quite common, Snakemake also has a more compact way of achieving the
same thing.

```python
input:
    "data/raw_internal/{sample_id}.fastq.gz",
    expand("intermediate/NCTC8325.{substr}.bt2",
        substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
```

> **Important!** <br>
> When using expand() like this, `substr` is not a wildcard because it is 
> resolved to the values explicitly given inside the expand expression.

Now change in the rules `index_genome` and `align_to_genome` to use the
`expand()` expression.

In the workflow we decide which samples to run by including the SRR ids in the
names of the inputs to the rules `multiqc` and `generate_count_table`:

```python
rule generate_count_table:
    output:
        "results/tables/counts.tsv"
    input:
        bams = ["intermediate/SRR935090.sorted.bam", 
                "intermediate/SRR935091.sorted.bam", 
                "intermediate/SRR935092.sorted.bam"],
...
rule multiqc:
    output:
        html = "results/multiqc.html",
        stats = "intermediate/multiqc_general_stats.txt"
    input:
        "intermediate/SRR935090_fastqc.zip",
        "intermediate/SRR935091_fastqc.zip",
        "intermediate/SRR935092_fastqc.zip"

```

The output files from these two rules, `results/multiqc.html` and  
`results/tables/counts.tsv`, are in turn specified as input to the `all` rule
at the top of the file. Because the first rule is targeted by default when 
we run Snakemake on the command line (like we mentioned in 
[snakemake-4-the-mrsa-workflow](snakemake-4-the-mrsa-workflow)) this 
is what triggers the rules to run on each of the three samples.

However, this is a potential source of errors since it's easy to change in one 
place and forget to change in the other. Because we can use Python code 
"everywhere" let's instead define a list of sample ids and put at the very 
top of the Snakefile, just before the rule `all`:

```python
SAMPLES = ["SRR935090", "SRR935091", "SRR935092"]
```

Now use `expand()` in `multiqc` and `generate_count_table` to use `SAMPLES` for
the sample ids. For the `multiqc` rule it could look like this:

```python
input:
    expand("intermediate/{sample_id}_fastqc.zip", sample_id = SAMPLES)
```

See if you can update the `generate_count_table` rule in the same manner!

> **Quick recap** <br>
> In this section we've learned:
>
> - How to use the `expand()` expression to create a list with file names, 
>   inserting all provided wildcard values.
