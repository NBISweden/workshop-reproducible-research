## MRSA workflow

This directory holds a final version of the MRSA Snakemake workflow structured
to be more in line with Snakemake best-practices.

It is meant as a **reference** that showcases how a Snakemake workflow can be
organised and written following the [Snakemake best
practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html).
It is _not_ the solution to the tutorial: the workflow you end up with after
going through the tutorial is simpler. To keep the scope of the tutorial
manageable, many of the features used here are not covered there (see [Features
not covered in the tutorial](#features-not-covered-in-the-tutorial) below).

### Folder structure

The workflow follows the [standard folder
structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility)
recommended for Snakemake workflows, where the workflow code lives under
`workflow/` and the project-specific settings live under `config/`:

```
final-workflow/
├── README.md                       # This file
├── environment.yml                 # Conda environment for running Snakemake itself
├── config/
│   ├── config.yml                  # Workflow settings (reads, genomes, sample sheet)
│   └── samples.csv                 # Sample sheet: sample ids and FASTQ URLs
└── workflow/
    ├── Snakefile                   # The workflow rules
    ├── envs/                       # One Conda environment file per tool
    │   ├── bowtie2.yml
    │   ├── fastqc.yml
    │   ├── multiqc.yml
    │   ├── samtools.yml
    │   ├── seqtk.yml
    │   └── subread.yml
    └── profiles/
        └── default/
            └── profile.yaml        # Default command line options for Snakemake
```

When the workflow runs it creates the following directories, which are not
tracked by git:

- `data/`: subsampled FASTQ files.
- `results/`: all outputs, including Bowtie2 indexes, BAM files, FastQC and
  MultiQC reports, the count table and log files under `results/logs/`.
- `.snakemake/`: Snakemake's metadata, the Conda environments and the files
  downloaded from remote storage (`.snakemake/storage/`).

### Configuration

- `config/config.yml` sets the number of reads to subsample per sample
  (`max_reads`), the path to the sample sheet (`sample_sheet`), the genome to
  align against (`genome_id`) and the URLs for the FASTA and GFF files of each
  available genome (`genomes`).
- `config/samples.csv` lists one sample per row, with the sample id and the URL
  to its FASTQ file.

To use another genome, add it under `genomes` and change `genome_id`. To use
other samples, update `samples.csv`.

### Workflow overview

1. `download_fastq_files` downloads the FASTQ file for each sample and
   subsamples it to `max_reads` reads with `seqtk`.
2. `fastqc` runs FastQC on each FASTQ file, and `multiqc` aggregates the FastQC
   reports into a single MultiQC report.
3. `index_genome` builds a Bowtie2 index of the genome selected by `genome_id`.
4. `align_to_genome` aligns each sample to the genome index, and `sort_bam`
   sorts the alignments into a BAM file.
5. `generate_count_table` counts the reads per gene for all samples with
   `featureCounts`.

The rule `all` collects the final outputs: `results/tables/counts.tsv` and
`results/multiqc/multiqc.html`.

### Running the workflow

Create and activate the Conda environment that provides Snakemake:

```bash
conda env create -f environment.yml
conda activate mrsa-workflow
```

Then run the workflow from this directory:

```bash
snakemake
```

Snakemake finds `workflow/Snakefile` automatically and loads the default
options in `workflow/profiles/default/profile.yaml`. This means one core, Conda
environments for each rule, the config file `config/config.yml`, keep going
when a job fails, and less verbose output. Any option given on the command line
overrides the profile, _e.g._ `snakemake -n` for a dry-run or `snakemake -c 4`
to use four cores. To use the containers defined in the rules instead of
Conda, run `snakemake --sdm apptainer` (this requires Apptainer).

### Features not covered in the tutorial

The following features are used in this workflow but are not, or only briefly,
covered in the tutorial:

- **Standard folder structure**: code in `workflow/`, settings in `config/`
  and one environment file per tool in `workflow/envs/`.
- **Workflow profile**: `workflow/profiles/default/profile.yaml` is loaded
  automatically and holds the default command line options, so the workflow
  runs with a plain `snakemake` command. The tutorial briefly covers profiles
  for running on HPC clusters.
- **Remote files via storage plugins**: the `storage:` directive and
  `storage.http()` let Snakemake download the FASTQ, FASTA and GFF files
  itself, before the jobs that use them run. This replaces the
  `curl`/`wget` commands and the `get_genome_fasta`/`get_genome_gff3` rules
  in the tutorial. See the [storage plugins in the plugin
  catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html).
- **Software environment for every rule**: each rule has a `conda:` directive
  pointing to an environment file with a pinned tool version, and rules
  also have a `container:` directive.
- **A log file for every rule**, including redirecting all output of a rule
  with `exec &>{log}`.
- **Named input and output files**, _e.g._ `{input.fastq}` and
  `{output.bam}` instead of `{input[0]}`.
- **The `multiext()` helper function** to define files that share a common
  prefix but have different extensions, here the Bowtie2 index files.
- **The `threads:` directive** to run Bowtie2 indexing with several threads.
- **The `message:` directive** to print a short description when a job
  starts.
- **`min_version()`** to require a minimum version of Snakemake.
- **Docstrings** that document the input functions.
- **A sample sheet read with pandas**. The tutorial covers this only in the
  extra material.

### Differences from the tutorial workflow

Apart from the features above, the workflow differs from the one built in the
tutorial in a few ways:

- There are no `get_genome_fasta` and `get_genome_gff3` rules. The genome files
  are downloaded through the storage plugin.
- Command line options, including which config file to use, are set in the
  profile instead of on the command line.
- Because each rule now has its dedicated conda environment file the main
  `environment.yml` file is slimmed down to just the packages required to run
  Snakemake. Any extra software requirements are handled by Snakemake.
