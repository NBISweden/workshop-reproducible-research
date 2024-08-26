We have up until now specified which Conda packages to install directly on the
command line using the `conda create` and `conda install` commands. For working
in projects this is not the recommended way. Instead, for increased control and
reproducibility, it is better to use an _environment file_ (in [YAML
format](https://en.wikipedia.org/wiki/YAML)) that specifies the packages,
versions and channels needed to create the environment for a project.

Throughout these tutorials we will use a case study where we analyse an RNA-seq
experiment with the multi-resistant bacteria MRSA (see [intro](introduction)).
You will now start to make a Conda YAML file for this MRSA project. The file
will contain a list of the software and versions needed to execute the analysis
code.

In this Conda tutorial, all code for the analysis is available in the script
`code/run_qc.sh`. This code will download the raw FASTQ-files and subsequently
run quality control on these using the FastQC software.

# Working with environments

We will start by making a Conda YAML-file that contains the required packages to
perform these two steps. Later in the course, you will update the Conda
YAML-file with more packages, as the analysis workflow is expanded.

- Let's get going! Make a YAML file called `environment.yml` looking like
  this, and save it in the current directory (which should be
  `workshop-reproducible-research/tutorials/conda`):

```yml
channels:
  - conda-forge
  - bioconda
dependencies:
  - fastqc=0.12.1
```

- Now, make a new Conda environment from the YAML file (note that here the
  command is `conda env create` as opposed to `conda create` that we used
  before):

```bash
conda env create -n project_mrsa -f environment.yml
```

> **Tip** <br>
> You can also specify exactly which channel a package should come from
> inside the environment file, using the `channel::package=version`
> syntax.

> **Tip** <br>
> Instead of the `-n` flag you can use the `-p` flag to set the full path to
> where the Conda environment should be installed. In that way you can
> contain the Conda environment inside the project directory, which does make
> sense from a reproducibility perspective, and makes it easier to keep track
> of what environment belongs to what project. If you don't specify `-p` the
> environment will be installed in the `envs/` directory inside your
> Conda installation path.

- Activate the environment!

- Now we can run the code for the MRSA project found in `code/run_qc.sh`,
  either by running `bash code/run_qc.sh` or by opening the `run_qc.sh` file
  and executing each line in the terminal one by one. Do this!

This should download the project FASTQ files and run FastQC on them (as
mentioned above).

- Check your directory contents (`ls -Rlh`, or in your file browser). It should
  now have the following structure:

```no-highlight
   conda/
    |
    |- code/
    |   |- run_qc.sh
    |
    |- data/
    |   |- SRR935090.fastq.gz
    |   |- SRR935091.fastq.gz
    |   |- SRR935092.fastq.gz
    |
    |- results/
    |   |- fastqc/
    |       |- SRR935090_fastqc.html
    |       |- SRR935090_fastqc.zip
    |       |- SRR935091_fastqc.html
    |       |- SRR935091_fastqc.zip
    |       |- SRR935092_fastqc.html
    |       |- SRR935092_fastqc.zip
    |
    |- environment.yml
```

Note that all that was needed to carry out the analysis and generate these
files and results was `environment.yml` (that we used to create a Conda
environment with the required packages) and the analysis code in
`code/run_qc.sh`.

# Keeping track of dependencies

Projects can often be quite large and require lots of dependencies; it can feel
daunting to try to capture all of that in a single Conda environment, especially
when you consider potential incompatibilities that may arise. It can therefore
be a good idea to start new projects with an environment file with each package
you know that you will need to use, but without specifying exact versions
(except for those packages where you _know_ you need a specific version).
This will install the latest compatible versions of all the specified
software, making the start-up and installation part of new projects easier. You
can then add the versions that were installed to your environment file
afterwards, ensuring future reproducibility.

There is one command that can make this easier: `conda env export`. This allows
you to export a list of the packages you've already installed, including their
specific versions, meaning you can easily add them after the fact to your
environment file. If you use the `--no-builds` flag, you'll get a list of the
packages minus their OS-specific build specifications, which is more useful for
making the environment portable across systems. This way, you can start with an
environment file with just the packages you need (without version), which will
install the most up-to-date version possible, and then add the resulting
version back in to the environment file using the `export` command!

> **Quick recap** <br>
> In this section we've learned:
>
> - How to define our Conda environment using a YAML-file.
> - How to use `conda env create` to make a new environment from a YAML-file.
> - How to use `conda env export` to get a list of installed packages.
> - How to work in a project-like setting.
