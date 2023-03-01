So far we have only defined the inputs/outputs of a rule as strings, or in
some case a list of strings, but Snakemake allows us to be much more flexible
than that. Actually, we can use any Python expression or even functions, as
long as they return a string or list of strings. Consider the rule
`align_to_genome` below.

```python
rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    output:
        "intermediate/{sample_id,\w+}.bam"
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

Here we have seven inputs; the fastq file with the reads and six files with
similar file names from the Bowtie 2 genome indexing. We can try to tidy this
up by using a Python expression to generate a list of these files instead. If
you're familiar with Python you could do this with
[list comprehensions](https://docs.python.org/3/tutorial/datastructures.html#list-comprehensions)
like this:

```python
input:
    fastq = "data/raw_internal/{sample_id}.fastq.gz",
    index = [f"intermediate/NCTC8325.{substr}.bt2" for
        substr in ["1", "2", "3", "4", "rev.1", "rev.2"]]
```

This will take the elements of the list of substrings one by one, and insert
that element in the place of `{substr}`. Since this type of aggregating
rules are quite common, Snakemake also has a more compact way of achieving the
same thing.

```python
input:
    fastq = "data/raw_internal/{sample_id}.fastq.gz",
    index = expand("intermediate/NCTC8325.{substr}.bt2",
           substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
```

Now change in the rules `index_genome` and `align_to_genome` to use the
`expand()` expression.

In the workflow we decide which samples to run by including the SRR ids in the
names of the inputs to the rules `multiqc` and `generate_count_table`. This is
a potential source of errors since it's easy to change in one place and forget
to change in the other. As we've mentioned before, but not really used so far,
Snakemake allows us to use Python code "everywhere". Let's therefore define
a list of sample ids and put at the very top of the Snakefile, just before the
rule `all`.

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

!!! Success "Quick recap"
    In this section we've learned:

    - How to use the `expand()` expression to create a list with file names,
    inserting all provided wildcard values.
