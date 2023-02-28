We have up until now specified which Conda packages to install directly on the
command line using the `conda create` and `conda install` commands. For working
in projects this is not the recommended way. Instead, for increased control and
reproducibility, it is better to use an *environment file* (in [YAML format](https://en.wikipedia.org/wiki/YAML))
that specifies the packages, versions and channels needed to create the
environment for a project.

Throughout these tutorials we will use a case study where we analyze an RNA-seq
experiment with the multiresistant bacteria MRSA (see [intro](introduction)).
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

* Let's get going! Make a YAML file called `environment.yml` looking like
  this, and save it in the current directory (which should be
  `workshop-reproducible-research/tutorials/conda`):

```yml
channels:
- conda-forge
- bioconda
dependencies:
- fastqc=0.11.9
- sra-tools=2.11.0
```

* Now, make a new Conda environment from the YAML file (note that here the
  command is `conda env create` as opposed to `conda create` that we used
  above):

```bash
conda env create -n project_mrsa -f environment.yml
```

!!! Tip
    You can also specify exactly which channel a package should come from
    inside the environment file, using the `channel::package=version`
    syntax.

!!! Tip
    Instead of the `-n` flag you can use the `-p` flag to set the full path to
    where the Conda environment should be installed. In that way you can
    contain the Conda environment inside the project directory, which does make
    sense from a reproducibility perspective, and makes it easier to keep track
    of what environment belongs to what project. If you don't specify `-p` the
    environment will be installed in the default `miniconda3/envs/` directory.

* Activate the environment!

* Now we can run the code for the MRSA project found in `code/run_qc.sh`,
  either by running `bash code/run_qc.sh` or by opening the `run_qc.sh` file
  and executing each line in the terminal one by one. Do this!

This should download the project FASTQ files and run FastQC on them (as
mentioned above).

* Check your directory contents (`ls -Rlh`, or in your file browser). It should
  now have the following structure:

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
(except for those packages where you *know* you need a specific version). Conda
will then try to get the latest compatible versions of all the specified
software, making the start-up and installation part of new projects easier. You
can then add the versions that were installed to your environment file
afterwards, ensuring future reproducibility.

There is one command that can make this easier: `conda env export`. This allows
you to export a list of the packages you've already installed, including their
specific versions, meaning you can easily add them after the fact to your
environment file. If you use the `--no-builds` flag, you'll get a list of the
packages minus their OS-specific build specifications, which is more useful for
making the environment portable across systems. This way, you can start with an
environment file with just the packages you need (without version), allow Conda
to solve the dependency tree and install the most up-to-date version possible,
and then add the resulting version back in to the environment file using the
`export` command!

# Optimising for speed

One of the greatest strengths of Conda is, unfortunately, also its greatest
weakness in its current implementation: the availability of a frankly enormous
number of packages and versions. This means that the search space for the
dependency hierarchy of any given Conda environment can become equally enormous,
leading to a (at times) ridiculous execution time for the dependency solver. It
is not uncommon to find yourself waiting for minutes for Conda to solve
a dependency hierarchy, sometimes even into the double digits. How can this be
circumvented?

Firstly, it is useful to specify as many of the `major.minor.patch` version
numbers as possible when defining your environment: this drastically reduces the
search space that Conda needs to go through. This is not always possible,
though. For example, we mentioned in the end of the *Environments in projects*
section that you might want to start out new projects without version
specifications for most packages, which means that the search space is going to
be large. Here is where another software comes into play: *Mamba*.

The [Mamba package manager](https://github.com/mamba-org/mamba) is built on-top
of Conda with some changes and additions that greatly speed up the execution
time. First of all, core parts of Mamba are written in C++ instead of Python,
like the original Conda. Secondly, it uses a different dependency solver
algorithm which is much faster than the one Conda uses. Lastly, it allows for
parallel downloading of repository data and package files with multi-threading.
All in all, these changes mean that Mamba is (currently) simply a better
version of Conda. Hopefully these changes will be incorporated into the Conda
core at some point in the future!

So, how do you get Mamba? Funnily enough, the easiest way to install it is (of
course) using Conda! Just run `conda install -n base -c conda-forge mamba`,
which will install Mamba in your `base` Conda environment. Mamba works almost
exactly the same as Conda, meaning that all you need to do is to stop using
`conda command` and instead use `mamba command` - simple! Be aware though that
in order to use `mamba activate` and `mamba deactivate` you first need to run
`mamba init`. So transitioning into using Mamba is actually quite easy - enjoy
your shorter execution times!

!!! Success "Quick recap"
    In this section we've learned:

    - How to define our Conda environment using a YAML-file.
    - How to use `conda env create` to make a new environment from a YAML-file.
    - How to use `conda env export` to get a list of installed packages.
    - How to work with Conda in a project-like setting.
    - How to optimise Conda for speed.
