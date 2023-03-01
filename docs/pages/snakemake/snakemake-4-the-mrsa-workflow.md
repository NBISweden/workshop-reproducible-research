As you might remember from the [intro](introduction), we are attempting to
understand how lytic bacteriophages can be used as a future therapy for the
multiresistant bacteria MRSA (methicillin-resistant _Staphylococcus aureus_).
In order to do this we have performed RNA-seq of three strains, one test and
two controls. We have already set up a draft Snakemake workflow for the RNA-seq
analysis and it seems to be running nicely. It's now up to you to modify this
workflow to make it more flexible and reproducible!

!!! Tip
    This section will leave a little more up to you compared to the previous
    one. If you get stuck at some point the final workflow after all the
    modifications is available in `tutorials/git/Snakefile`.

You are probably already in your `snakemake-env` environment, otherwise
activate it (use `conda info --envs` if you are unsure).

!!! Tip
    Here we have one Conda environment for executing the whole Snakemake
    workflow. Snakemake also supports using explicit Conda environments on
    a per-rule basis, by specifying something like `conda:
    rule-specific-env.yml` in the rule definition and running Snakemake with
    the `--use-conda` flag. The given rule will then be run in the Conda
    environment specified in `rule-specific-env.yml` that will be created and
    activated on the fly by Snakemake.

Let's start by generating the rule graph so that we get an overview of the
workflow.

```bash
snakemake -s snakefile_mrsa.smk --rulegraph | dot -T png > rulegraph_mrsa.png
```

There are two differences in this command compared to the one we've used
before. The first is that we're using the `-s` flag to specify which Snakemake
workflow to run. We didn't need to do that before since `Snakefile` is the
default name. The second is that we don't define a target. In the toy example
we used `a_b.txt` as a target, and the wildcards were resolved based on that.
How come that we don't need to do that here? It turns out that by default
Snakemake targets the first rule in a workflow. By convention, we call this rule
`all` and let it serve as a rule for aggregating the main outputs of the
workflow.

![](images/rulegraph_mrsa.svg)

Now take some time and look through the workflow file and try to understand how
the rules fit together. Use the rule graph as aid. The rules represent a quite
standard, although somewhat simplified, workflow for RNA-seq analysis. If you
are unfamiliar with the purpose of the different operations (index genome,
FastQC and so on), then take a look at the [intro](introduction).

Also generate the job graph in the same manner. Here you can see that three
samples will be downloaded from SRA (Sequence Read Archive); SRR935090,
SRR935091, and SRR935092. Those will then be quality controlled with FastQC and
aligned to a genome. The QC output will be aggregated with MultiQC and the
alignments will be used to generate a count table, *i.e.* a table that shows
how many reads map to each gene for each sample. This count table is then what
the downstream analysis will be based on.

![](images/dag_mrsa.svg)

Now try to run the whole workflow. Hopefully you see something like this.

```no-highlight
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job stats:
job                     count    min threads    max threads
--------------------  -------  -------------  -------------
align_to_genome             3              1              1
all                         1              1              1
fastqc                      3              1              1
generate_count_table        1              1              1
generate_rulegraph          1              1              1
get_SRA_by_accession        3              1              1
get_genome_fasta            1              1              1
get_genome_gff3             1              1              1
index_genome                1              1              1
multiqc                     1              1              1
sort_bam                    3              1              1
total                      19              1              1

Select jobs to execute...

[Mon Oct 25 17:13:47 2021]
rule get_genome_fasta:
    output: data/raw_external/NCTC8325.fa.gz
    jobid: 6
    resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T

--2021-10-25 17:13:48--  ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
           => ‘data/raw_external/NCTC8325.fa.gz’
Resolving ftp.ensemblgenomes.org (ftp.ensemblgenomes.org)... 193.62.197.75
Connecting to ftp.ensemblgenomes.org (ftp.ensemblgenomes.org)|193.62.197.75|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
.
.
[lots of stuff]
.
.
localrule all:
    input: results/tables/counts.tsv, results/multiqc.html, results/rulegraph.png
    jobid: 0
    resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T

[Mon Oct 25 17:14:38 2021]
Finished job 0.
19 of 19 steps (100%) done
```

After everything is done, the workflow will have resulted in a bunch of files
in the directories `data`, `intermediate` and `results`. Take some time to look
through the structure, in particular the quality control reports in `results`
and the count table in `results/tables`.

!!! Success "Quick recap"
    In this section we've learned:

    - How the MRSA workflow looks.
    - How to run the MRSA workflow.
    - Which output files the MRSA workflow produces.
