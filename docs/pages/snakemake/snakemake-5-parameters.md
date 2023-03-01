In a typical bioinformatics project, considerable efforts are spent on tweaking
parameters for the various programs involved. It would be inconvenient if you
had to change in the shell scripts themselves every time you wanted to run with
a new setting. Luckily, there is a better option for this: the `params`
keyword.

```python
rule some_rule:
    output:
        "..."
    input:
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
(`-X 25000`). Change in this rule to use the parameter `max_reads` instead and
set the value to 20000. If you need help, click to show the solution below.

??? example "Click to show the solution"
    ```python
    rule get_SRA_by_accession:
        """
        Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
        """
        output:
            "data/raw_internal/{sample_id}.fastq.gz"
        params:
            max_reads = 20000
        shell:
            """
            fastq-dump {wildcards.sample_id} -X {params.max_reads} --readids \
                --dumpbase --skip-technical --gzip -Z > {output}
            """
    ```

Now run through the workflow. Because there's been changes to the `get_SRA_by_accession`
rule this will trigger a re-run of the rule for all three accessions. In addition
all downstream rules that depend on output from `get_SRA_by_accession` are re-run.

The parameter values we set in the `params` section don't have to be static,
they can be any Python expression. In particular, Snakemake provides a global
dictionary of configuration parameters called `config`. Let's modify
`get_SRA_by_accession` to look something like this in order to make use of this
dictionary:

```python
rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sample_id}.fastq.gz"
    params:
        max_reads = config["max_reads"]
    shell:
        """
        fastq-dump {wildcards.sample_id} -X {params.max_reads} --readids \
            --dumpbase --skip-technical --gzip -Z > {output}
        """
```

Note that Snakemake now expects there to be a key named `max_reads` in the config
dictionary. If we don't populate the dictionary somehow the dictionary will be
empty so if you were to run the workflow now it would trigger a `KeyError` (try
running `snakemake -s snakefile_mrsa.smk -n` to see for yourself).
In order to populate the config dictionary with data for the workflow we could
use the `snakemake --config KEY=VALUE` syntax directly from the command line.
However, from a reproducibility perspective, it's not optimal to set parameters
from the command line, since it's difficult to keep track of which parameter
values that were used.

A much better alternative is to use the `--configfile FILE` option to supply a
configuration file to Snakemake. In this file we can collect all the
project-specific settings, sample ids and so on. This also enables us to write
the Snakefile in a more general manner so that it can be better reused between
projects. Like several other files used in these tutorials, this file should be
in [yaml format](https://en.wikipedia.org/wiki/YAML). Create the file below and
save it as `config.yml`.

```yaml
max_reads: 25000
```

If we now run Snakemake with `--configfile config.yml`, it will parse this file
to form the `config` dictionary. If you want to overwrite a parameter value,
*e.g.* for testing, you can still use the `--config KEY=VALUE` flag, as in
`--config max_reads=1000`.

!!! Tip
    Rather than supplying the config file from the command line you could also
    add the line `configfile: "config.yml"` to the top of your Snakefile. Keep in
    mind that with such a setup Snakemake will complain if the file `config.yml`
    is not present.

!!! Success "Quick recap"
    In this section we've learned:

    - How to set parameter values with the `params` directive.
    - How to run Snakemake with the `config` variable and with a configuration file.
