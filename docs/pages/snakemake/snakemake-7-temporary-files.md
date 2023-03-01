It's not uncommon that workflows contain temporary files that should be kept
for some time and then deleted once they are no longer needed. A typical case
could be that some operation generates a file, which is then compressed to save
space or indexed to make searching faster. There is then no need to save the
original output file. Take a look at the job graph for our workflow again. The
output from `align_to_genome` is a bam file, which contains information about
all the reads for a sample and where they map in the genome. For downstream
processing we need this file to be sorted by genome coordinates. This is what
the rule `sort_bam` is for. We therefore end up with both
`intermediate/{sample_id}.bam` and `intermediate/{sample_id}.sorted.bam`.

In Snakemake we can mark an output file as temporary like this:

```python
output: temp("...")
```

The file will then be deleted as soon as all jobs where it's an input have
finished. Now do this for the output of `align_to_genome`. We have to rerun the
rule for it to trigger, so use `-R align_to_genome`. It should look something
like this:

```no-highlight
.
.
rule sort_bam:
    input: intermediate/SRR935090.bam
    output: intermediate/SRR935090.sorted.bam
    jobid: 2
    wildcards: sample_id=SRR935090

Removing temporary output file intermediate/SRR935090.bam.
Finished job 2.
.
.
```

!!! Tip
    Sometimes you may want to trigger removal of temporary files without
    actually rerunning the jobs. You can then use the `--delete-temp-output`
    flag. In some cases you may instead want to run only parts of a workflow
    and therefore want to prevent files marked as temporary from being deleted
    (because the files are needed for other parts of the workflow). In such
    cases you can use the `--notemp` flag.

Snakemake has a number of options for marking files:

* `temp("...")`: The output file should be deleted once it's no longer needed
  by any rules.
* `protected("...")`: The output file should be write-protected. Typically used
  to protect files that require a huge amount of computational resources from
  being accidentally deleted.
* `ancient("...")`: The timestamp of the input file is ignored and it's always
  assumed to be older than any of the output files.
* `touch("...")`: The output file should be "touched", *i.e.* created or
  updated, when the rule has finished. Typically used as "flag files" to
  enforce some rule execution order without real file dependencies.
* `directory("...")`: The output is a directory rather than a file.

Success "Quick recap"
    In this section we've learned:

    - How to mark an output file as temporary for automatic removal.
