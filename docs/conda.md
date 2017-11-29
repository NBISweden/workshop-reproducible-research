# Conda tutorial

## Introduction
### What is Conda?
Conda is a package and environment manager. As a package manager it enables you to install a wide range of software and tools using one simple command (`conda install`). As an environment manager it allows you to create and manage multiple different environments, each with their own set of packages. Why would you want to do that? For instance, you may want to be able to easily run different versions of the same package or have different cross-package dependencies that are incompatible with each other.

Environments are of particular relevance when making bioinformatics projects reproducible. Full reproducibility requires the possibility to recreate the system that was originally used to generate the results. This can (to a large extent) be accomplished by using Conda to make a project environment with specific versions of a set of packages that are needed to execute and run all code needed in the project. You can read more about Conda [here](https://conda.io/docs/user-guide/concepts.html).

### Conda packages
A Conda package is a compressed tarball (system-level libraries, Python or other modules, executable programs, or other components). Conda keeps track of the dependencies between packages and platforms - this means that when installing a given package, all necessary dependencies will also be installed. Conda packages are hosted and downloaded from remote channels. Some widely used channels for general-purpose and bioinformatics packages are [Condaforge](https://conda-forge.org/) and [Bioconda](https://bioconda.github.io/). When installing a Conda package you specify the package name, version (optional), and channel to download from.

### Conda environments
A Conda environment is essentially a directory that contains a specific collection of conda packages that you have installed. Packages are symlinked between environments to avoid unnecessary duplication.

## Practical exercise

This tutorial depends on files from the course BitBucket repo. Take a look at the [intro](index) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/conda`. Instructions below assume that you are standing in `conda/` unless otherwise specified (e.g. if it says "create a file", it means save it in `conda/`).

### Install Conda (not for Windows)
If you are doing the tutorials by running a Docker container on your Windows machine, Conda will already be installed for you. You can jump down to the next section (Get going with environments).


The easiest way to get going is to install Conda by downloading the correct binary from [here](https://conda.io/miniconda.html) and running the installer. This can be done on the command line as shown below:
```bash
# --- Example (don't run) ---
# Install Miniconda on 64-bit Linux for python3:
wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh
```

* Go ahead and install Conda as described above, but be sure to modify the `wget` command to download the correct file for your system (see link above).
* The installer will ask you questions during the installation:
  - do you accept the license terms?
  - do you accept the installation path or do you want to chose a different one?
  - do you want to add the Conda installation to your PATH by adding a line to your `.bashrc` file?  
   *If you select yes, the system will find Conda and you can run the `conda` command directly. If you select no, you will either have to run Conda by specifying the full path, e.g. `/root/miniconda3/bin/conda`, or by adding it manually to your PATH each time you restart your system, e.g. `export PATH="/root/miniconda3/bin:$PATH"`.*
* You can verify the version you installed by running:

```bash
conda --version
```

### Get going with environments
* Let's make our first Conda environment:

```bash
conda create -n Project_A -c bioconda fastqc
```
This will create an environment called Project_A containing FastQC from the Bioconda channel. Conda will list the packages that will be installed and ask for your confirmation.
* Once it is done, you can activate the environment:

```
source activate Project_A
```
By default, Conda will add information to your prompt telling you which environment that is active.
* To see all your environments you can run:

```
conda info --envs
```
The active environment will be marked with an asterisk.
* To see the installed packages in the active environment, run:

```
conda list
```
* Now, deactivate the environment by running `source deactivate`.
* List all environments again. Which environment is now marked as active?
* Try to run FastQC:

```
fastqc --version
```
* Did it work? Activate your Project_A environment and run the `fastqc --version` command again. Does it work now?

Hopefully the FastQC software was not found in your root environment (unless you had installed it previously) but worked once your environment was activated!
* Now, let's add another package (sra-tools) to our environment using `conda install`. Make sure that Project_A is the active environment first.

```
conda install -c bioconda sra-tools
```
* If we don't specify the package version, the latest available version will be installed. What version of sra-tools got installed?
* Run the following to see what versions are available:

```
conda search -c bioconda sra-tools
```

* Now try to install a different version of sra-tools, e.g.:

```
conda install -c bioconda sra-tools=2.7.0
```

Read the information that Conda displays in the terminal. It probably asks if you want to downgrade the initial sra-tools installation to the one specified here (2.7.0 in the example). You can only have one version of a given package in a given environment.

Let's assume that you are just about to start a new exciting research project called Project A. You will have sequencing data and want to use the latest bowtie2 software to align your reads.
* Find out what versions of bowtie2 are available at Bioconda using `conda search`. Now install the *latest* available version of bowtie2 in your Project_A environment.

Let's further assume that you have an old project (called Project Old) where you know you used bowtie2 v2.2.6. You just got back reviewer comments and they want you to include some alignment statistics. Unfortunately you haven't saved that information so you will have to rerun the alignment. Now, it is essential that you use the same version of bowtie that your results are based on, otherwise the alignment statistics will be misleading. Using Conda environments this becomes simple. You can just have a separate environment for your old project where you have an old version of bowtie2 without interfering with your new Project A where you want the latest version.

* Make a new environment for your old project:

```
conda create -n Project_Old -c bioconda bowtie2=2.2.6
```
* Activate `Project_Old` and check the bowtie2 version (`bowtie2 --version`).
* Activate `Project_A` again and check the bowtie2 version.
* Run `source deactivate` to exit your active environment.
* List your environments (do you remember the command?).
* Now, let's remove an environment:

```
conda env remove -n Project_Old
```
After making a few different environments and installing a bunch of packages, Conda can take up some disk space. In order to remove unnecessary files one can use the command:

```
conda clean -a
```
This will remove package tar-balls that are left from package installations, unused packages (i.e. not present in any environments), and cached data.

### How to use in a reproducible project setting

We have up until now specified which Conda packages to install directly on the command line using the `conda create` and `conda install` commands. For working in projects this is not the recommended way. Instead, for increased control and reproducibility, it is better to use a file (in [yaml format](https://en.wikipedia.org/wiki/YAML)) specifying packages, versions and channels needed to create the environment for a project.

* You will now start to make a Conda yml-file for the MRSA project that we will use in this course!
* To start with, the analysis code is contained in one file, `code/run_qc.sh`, and will:
  - download the raw fastq-files
  - run quality control on these using FastQC
* We will start by making a Conda yml-file that contains the required packages to perform these steps.
* Later in the course, you will update the Conda yml-file with more packages, as the analysis workflow is expanded.
* Make a yaml file called `environment.yml` for our tutorial project looking like this:

```yaml
channels:
- conda-forge
- bioconda
dependencies:
- fastqc=0.11
- sra-tools=2.8
```

* Now, make a new Conda environment from the yaml file:

```bash
conda env create -n Project_MRSA -f environment.yml
```
* Activate the environment!
* Now we can run the MRSA project analysis code found in `code/run_qc.sh`, either by

```bash
bash code/run_qc.sh
```
or by opening the `run_qc.sh` file and executing each line in the terminal one by one. Do this!

This should download the project fastq files and run FastQC on them (as mentioned above).
* Check your directory content (`ls -Rlh`, or in your file browser). It should now look something like this:

```
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
Note that all that was needed to carry out the analysis and generate these files and results was the `environment.yml` file (that we used to create a Conda environment with the required packages) and the analysis code in `code/run_qc.sh`.

## Where to go next?

### Version tracking with git
It is obvious that these two files (`environment.yml` and `run_qc.sh`) are important to not lose as they are essential to reproduce the analysis. A very convenient way to keep track of these files (and others) is to set up git version control. The content of these files, and their history, can then be safely tracked and stored in the cloud (e.g. at GitHub or Bitbucket). [Go to the git tutorial](git).

### Workflow management
We now performed several analysis steps through the code in `run_qc.sh`. That is ok, but it can quickly get complicated when the analysis workflow gets more complex and is split into several code files. Suddenly it is not obvious which code to run and in what order it should be executed, in order to rerun the full project analysis from scratch. A good solution to this is to use a workflow manager. [Go to the Snakemake tutorial](snakemake).
