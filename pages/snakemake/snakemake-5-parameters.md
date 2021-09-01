In a typical bioinformatics project, considerable efforts are spent on tweaking
parameters for the various programs involved. It would be inconvenient if you
had to change in the shell scripts themselves every time you wanted to run with
a new setting. Luckily, there is a better option for this: the `params`
keyword.

```python
rule some_rule:
    input:
        "..."
    output:
        "..."
    params:
        cutoff=2.5
    shell:
        """
        some_program --cutoff {params.cutoff} {input} {output}
        """
```

We run most of the programs with default settings in our workflow. However,
there is one parameter in the rule `get_SRA_by_accession` that we use for
determining how many reads we want to retrieve from SRA for each sample 
(`-X 25000`). Change in this rule to use the parameter `max_reads` instead, set the
value to 20000, and run through the workflow. Remember that Snakemake doesn't
automatically rerun rules after parameter changes, so you have to trigger the
execution of `get_SRA_by_accession` with `-R`.

```python
rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sra_id}.fastq.gz"
    shell:
        """
        fastq-dump {wildcards.sra_id} -X 25000 --readids \
            --dumpbase --skip-technical --gzip -Z > {output}
        """
```

The parameter values we set in the `params` section don't have to be static,
they can be any Python expression. In particular, Snakemake provides a global
dictionary of configuration parameters called `config`. Let's modify
`get_SRA_by_accession` to look something like this in order to access the
elements of this dictionary:

```python
rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sra_id}.fastq.gz"
    params:
        max_reads = config["max_reads"]
    shell:
        """
        fastq-dump {wildcards.sra_id} -X {params.max_reads} --readids \
            --dumpbase --skip-technical --gzip -Z > {output}
        """
```

The `config` variable is just a normal Python dictionary, but it has the
special feature that we can change the parameter values from the command line
by using the `snakemake --config KEY=VALUE` syntax. Try this out for yourself.

From a reproducibility perspective, it's not optimal to set parameters from the
command line, since it's difficult to keep track of which parameter values that
were used. A much better alternative is to use the `--configfile FILE` option.
Here we can collect all the project-specific settings, sample ids, and so on in
one file. This also enables us to write the Snakefile in a more general manner
so that it can be better reused between projects. Like several other files used
in these tutorials, this file should be in [yaml
format](https://en.wikipedia.org/wiki/YAML). Create the file below and save it
as `config.yml`.

```yaml
max_reads: 25000
```

If we now run Snakemake with `--configfile config.yml`, it will parse this file
to form the `config` dictionary. If you want to overwrite a parameter value,
*e.g.* for testing, you can still use the `--config` flag.

> **Tip** <br>
> Rather than supplying the config file from the command line you could also
> add the line `configfile: "config.yml"` to the top of your Snakefile.
