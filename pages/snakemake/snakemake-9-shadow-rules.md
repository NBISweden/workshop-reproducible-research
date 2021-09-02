Take a look at the rule `generate_count_table` below. Since `input.annotation`
is compressed, it is first unzipped to a temporary file. `htseq-count` then
generates a temporary count table, which is finally prepended with a header
containing the sample names. Lastly, the two temporary files are deleted.

```python
rule generate_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams=["intermediate/SRR935090.sorted.bam", "intermediate/SRR935091.sorted.bam", "intermediate/SRR935092.sorted.bam"],
        annotation="data/raw_external/NCTC8325.gff3.gz"
    output:
        "results/tables/counts.tsv"
    shell:
        """
        # htseq-count cannot use .gz, so unzip to a temporary file first
        gunzip -c {input.annotation} > tempfile

        # Save the count table as a temporary file and then prepend a header line
        # with the sample names
        htseq-count --format bam --type gene --idattr gene_id {input.bams} tempfile > tempfile2
        echo '{input.bams}' | tr ' ' '\t' | cat - tempfile2 > {output}

        # Remove the temporary files
        rm tempfile tempfile2
        """
```

There are a number of drawbacks with having files that aren't explicitly part
of the workflow as input/output files to rules (as the two temporary files
here).

* Snakemake cannot clean up these files if the job fails, as it would do for
  normal output files.
* If several jobs are run in parallel there is a risk that they write to
  `tempfile` at the same time. This can lead to very scary results.
* Sometimes we don't know the names of all the files that a program can
  generate. It is, for example, not unusual that programs leave some kind of
  error log behind if something goes wrong.

All of these issues can be dealt with by using the `shadow` option for a rule.
The shadow option results in that each execution of the rule is run in an
isolated temporary directory (located in `.snakemake/shadow/` by default).
There are a few options for `shadow`. The most simple is `shadow: "minimal"`,
which means that the rule is executed in an empty directory that the input
files to the rule have been symlinked into. For the rule below, that means that
the only file available would be `input.txt`. The shell commands would generate
the files `some_other_junk_file` and `output.txt`. Lastly, Snakemake will move
the output file (`output.txt`) to its "real" location and remove the whole
shadow directory. We therefore never have to think about manually removing
`some_other_junk_file`.

```python
rule some_rule:
    input:
        "input.txt"
    output:
        "output.txt"
    shadow: "minimal"
    shell:
        """
        touch some_other_junk_file
        cp {input} {output}
        """
```

Try this out for the rules where we have to "manually" deal with files that
aren't tracked by Snakemake (`multiqc`, `index_genome`,
`generate_count_table`). Also remove the shell commands that remove temporary
files from those rules, as they are no longer needed. Now rerun the workflow
and validate that the temporary files don't show up in your working directory.

> **Tip** <br>
> Some people use the shadow option for almost every rule and some never
> use it at all. One thing to keep in mind is that it leads to some extra file
> operations when the outputs are moved to their final location. This is no
> issue when the shadow directory is on the same disk as the output directory,
> but if you're running on a distributed file system and generate very many
> or very large files it might be worth considering other options (see *e.g.*
> the `--shadow-prefix` flag).

> **Quick recap** <br>
> In this section we've learned:
>
> - How to use the shadow option to handle files that are not tracked by Snakemake.