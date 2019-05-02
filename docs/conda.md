# Introduction to Conda
## What is Conda?
Conda is a package and environment manager. As a package manager it enables you to install a wide range of software and tools using one simple command: `conda install`. As an environment manager it allows you to create and manage multiple different environments, each with their own set of packages. Why would you want to do that? For instance, you can easily run different versions of the same package or have different cross-package dependencies that are otherwise incompatible with each other.

Environments are of particular relevance when making bioinformatics projects reproducible. Full reproducibility requires the possibility to recreate the system that was originally used to generate the results. This can, to a large extent, be accomplished by using Conda to make a project environment with specific versions of the packages that are needed in the project. You can read more about Conda [here](https://conda.io/docs/user-guide/concepts.html).

## Conda packages
A Conda package is a compressed tarball (system-level libraries, Python or other modules, executable programs, or other components). Conda keeps track of the dependencies between packages and platforms - this means that when installing a given package, all necessary dependencies will also be installed. Conda packages are typically hosted and downloaded from remote channels. Some widely used channels for general-purpose and bioinformatics packages are [conda-forge](https://conda-forge.org/) and [Bioconda](https://bioconda.github.io/), respectively. Both of these are community-driven projects, so if you're missing some package you can contribute by adding it to either channel. When installing a Conda package you specify the package name, version (optional), and channel to download from.

## Conda environments
A Conda environment is essentially a directory that is added to your PATH and that contains a specific collection of Conda packages that you have installed. Packages are symlinked between environments to avoid unnecessary duplication.

# Set up
This tutorial depends on files from the course Bitbucket repo. Take a look at the [intro](tutorial_intro.md) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/conda`. Instructions below assume that you are standing in `conda/` unless otherwise specified (e.g. if it says "create a file", it means save it in `conda/`).

## Install Conda
!!! attention
    If you are doing the tutorials by running a Docker container on your Windows machine, Conda will already be installed for you. You can then skip this section and go directly to the practical exercise.

!!! attention
    If you already have installed conda but want to update, you should be able to simply run `conda update conda` and subsequently `conda init`, and skip the installation instructions below.

* Go ahead and install Conda as described below. **Be sure to download the correct file for your OS**.

```bash
# *** Install Miniconda3 for 64-bit macOS ***
curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-MacOSX-x86_64.sh -O
bash Miniconda3-4.6.14-MacOSX-x86_64.sh
rm Miniconda3-4.6.14-MacOSX-x86_64.sh

# *** Install Miniconda3 for 64-bit Linux ***
curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O
bash Miniconda3-4.6.14-Linux-x86_64.sh
rm Miniconda3-4.6.14-Linux-x86_64.sh
```

* The installer will ask you questions during the installation:
    - do you accept the license terms? (Yes)
    - do you accept the installation path or do you want to chose a different one? (Probably yes)
    - do you want to run conda init to setup conda on your system? (Yes)
* Either restart your shell so the settings in `~/.bashrc`/`~/.bash_profile` can take effect, or source `~/.bashrc`/`~/.bash_profile`.
* You can verify that the installation worked by running:

```bash
conda --version
```

* Next, we will setup the the default channels (from where packages will be searched for and downloaded if no channel is specified).

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

# Practical exercise

## Get going with environments

Let's assume that you are just about to start a new exciting research project called Project A.

* Let's make our first Conda environment:

```bash
conda create -n project_a -c bioconda fastqc
```

This will create an environment called project_a, containing FastQC from the Bioconda channel. Conda will list the packages that will be installed and ask for your confirmation.

* Once it is done, you can activate the environment:

```bash
conda activate project_a
```

By default, Conda will add information to your prompt telling you which environment that is active.

* To see all your environments you can run:

```bash
conda info --envs
```

The active environment will be marked with an asterisk.

* To see the installed packages in the active environment, run:

```bash
conda list
```

* Now, deactivate the environment by running `conda deactivate`.
* List all environments again. Which environment is now marked as active?
* Try to run FastQC:

```bash
fastqc --version
```

* Did it work? Activate your project_a environment and run the `fastqc --version` command again. Does it work now?

Hopefully the FastQC software was not found in your base environment (unless you had installed it previously), but worked once your environment was activated.

* Now, let's add another package (sra-tools) to our environment using `conda install`. Make sure that project_a is the active environment first.

```bash
conda install -c bioconda sra-tools
```

* If we don't specify the package version, the latest available version will be installed. What version of sra-tools got installed?
* Run the following to see what versions are available:

```bash
conda search -c bioconda sra-tools
```

* Now try to install a different version of sra-tools, e.g.:

```bash
conda install -c bioconda sra-tools=2.7.0
```

Read the information that Conda displays in the terminal. It probably asks if you want to downgrade the initial sra-tools installation to the one specified here (2.7.0 in the example). You can only have one version of a given package in a given environment.

Let's assume that you will have sequencing data in your Project A, and want to use the latest bowtie2 software to align your reads.

* Find out what versions of bowtie2 are available in the bioconda channel using `conda search -c bioconda`.
* Now install the *latest* available version of bowtie2 in your project_a environment.

Let's further assume that you have an old project (called Project Old) where you know you used bowtie2 v2.2.5. You just got back reviewer comments and they want you to include some alignment statistics. Unfortunately, you haven't saved that information so you will have to rerun the alignment. Now, it is essential that you use the same version of bowtie that your results are based on, otherwise the alignment statistics will be misleading. Using Conda environments this becomes simple. You can just have a separate environment for your old project where you have an old version of bowtie2 without interfering with your new Project A where you want the latest version.

* Make a new environment for your old project:

```bash
conda create -n project_old -c bioconda bowtie2=2.2.5
```

* Activate `project_old` and check the bowtie2 version (`bowtie2 --version`).
* Activate `project_a` again and check the bowtie2 version.
* Run `conda deactivate` to exit your active environment.
* List your environments (do you remember the command?).
* Now, let's remove an environment:

```bash
conda env remove -n project_old
```

After making a few different environments and installing a bunch of packages, Conda can take up some disk space. You can remove unnecessary files with the command:

```bash
conda clean -a
```

This will remove package tar-balls that are left from package installations, unused packages (i.e. those not present in any environments), and cached data.

!!! note "Quick recap"
    In this section we've learned:

    * How to use `conda install` for installing packages.
    * How to create and activate environments and how to change between them.
    * How to remove packages or environments and clean up.

## How to use in a reproducible project setting

We have up until now specified which Conda packages to install directly on the command line using the `conda create` and `conda install` commands. For working in projects this is not the recommended way. Instead, for increased control and reproducibility, it is better to use a file (in [yaml format](https://en.wikipedia.org/wiki/YAML)) specifying packages, versions and channels needed to create the environment for a project.

* Throughout these tutorials we will use a case study where we analyze an RNA-seq experiment with the multiresistant bacteria MRSA (see [intro](tutorial_intro.md)). You will now start to make a Conda yml file for this MRSA project. The file will contain a list of the software and versions needed to execute the analysis code.
* In this conda tutorial, all code for the analysis is available in the script `code/run_qc.sh`. This code will:
    - download the raw fastq-files
    - run quality control on these using FastQC
* We will start by making a Conda yml-file that contains the required packages to perform these two steps.
* Later in the course, you will update the Conda yml-file with more packages, as the analysis workflow is expanded.
* So let's get going! Make a yaml file called `environment.yml` looking like this, and save it in the current directory (should be `reproducible_research_course/conda`):

```yaml
channels:
- conda-forge
- bioconda
dependencies:
- fastqc=0.11.6
- sra-tools=2.8
```

* Now, make a new Conda environment from the yaml file:

```bash
conda env create -n project_mrsa -f environment.yml
```

* Activate the environment!
* Now we can run the code for the MRSA project found in `code/run_qc.sh`, either by running `bash code/run_qc.sh` or by opening the `run_qc.sh` file and executing each line in the terminal one by one. Do this!

This should download the project fastq files and run FastQC on them (as mentioned above).

* Check your directory contents (`ls -Rlh`, or in your file browser). It should now have the following structure:

```no-highlight
   conda/
    |
    |- code/
    |   |- run_qc.sh
    |
    |- data/
    |   |- raw_internal/
    |       |- SRR935090.fastq.gz
    |       |- SRR935091.fastq.gz
    |       |- SRR935092.fastq.gz
    |
    |- intermediate/
    |   |- fastqc/
    |       |- SRR935090_fastqc.zip
    |       |- SRR935091_fastqc.zip
    |       |- SRR935092_fastqc.zip
    |
    |- results/
    |   |- fastqc/
    |       |- SRR935090_fastqc.html
    |       |- SRR935091_fastqc.html
    |       |- SRR935092_fastqc.html
    |
    |- environment.yml
```

Note that all that was needed to carry out the analysis and generate these files and results was `environment.yml` (that we used to create a Conda environment with the required packages) and the analysis code in `code/run_qc.sh`.

!!! note "Quick recap"
    In this section we've learned:

    * How to define our Conda environment using a yaml-file.
    * How to use `conda env create` to make a new environment from a yaml-file.
