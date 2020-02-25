# Introduction to Conda

## What is Conda?
Conda is a package and environment manager. As a package manager it enables you
to install a wide range of software and tools using one simple command: `conda
install`. As an environment manager it allows you to create and manage multiple
different environments, each with their own set of packages. Why would you want
to do that? For instance, you can easily run different versions of the same
package or have different cross-package dependencies that are otherwise
incompatible with each other.

Environments are of particular relevance when making bioinformatics projects
reproducible. Full reproducibility requires the possibility to recreate the
system that was originally used to generate the results. This can, to a large
extent, be accomplished by using Conda to make a project environment with
specific versions of the packages that are needed in the project. You can read
more about Conda [here](https://conda.io/docs/user-guide/concepts.html).

## Conda packages
A Conda package is a compressed tarball (system-level libraries, Python or
other modules, executable programs, or other components). Conda keeps track of
the dependencies between packages and platforms - this means that when
installing a given package, all necessary dependencies will also be installed.
Conda packages are typically hosted and downloaded from remote channels. Some
widely used channels for general-purpose and bioinformatics packages are
[conda-forge](https://conda-forge.org/) and
[Bioconda](https://bioconda.github.io/), respectively. Both of these are
community-driven projects, so if you're missing some package you can contribute
by adding it to either channel. When installing a Conda package you specify the
package name, version (optional), and channel to download from.

## Conda environments
A Conda environment is essentially a directory that is added to your PATH and
that contains a specific collection of Conda packages that you have installed.
Packages are symlinked between environments to avoid unnecessary duplication.

# Set up
This tutorial depends on files from the course GitHub repo. Take a look at the
[intro](tutorial_intro.md) for instructions on how to set it up if you haven't
done so already. Then open up a terminal and go to
`workshop-reproducible-research/conda`. Instructions below assume that you are
standing in `conda/` unless otherwise specified (*e.g.* if it says "create
a file", it means save it in `conda/`).

## Install Conda
!!! attention
    If you are doing the tutorials by running a Docker container on your
    Windows machine, Conda will already be installed for you. You can then jump
    ahead to the last point about setting up the default channels (`conda
    config`) and then go ahead with the practical exercise.

!!! attention
    If you already have installed conda but want to update, you should be able
    to simply run `conda update conda` and subsequently `conda init`, and skip
    the installation instructions below.

* Go ahead and install Conda as described below. **Be sure to download the
  correct file for your OS**.

**macOS**
```bash
# *** Install Miniconda3 for 64-bit macOS ***
curl https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-MacOSX-x86_64.sh -O
bash Miniconda3-4.7.12.1-MacOSX-x86_64.sh
rm Miniconda3-4.7.12.1-MacOSX-x86_64.sh
```

**Linux**
```bash
# *** Install Miniconda3 for 64-bit Linux ***
curl https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O
bash Miniconda3-4.7.12.1-Linux-x86_64.sh
rm Miniconda3-4.7.12.1-Linux-x86_64.sh
```

* The installer will ask you questions during the installation:
    - Do you accept the license terms? (Yes)
    - Do you accept the installation path or do you want to chose a
        different one? (Probably yes)
    - Do you want to run conda init to setup conda on your system? (Yes)
* Restart your shell so the settings in `~/.bashrc`/`~/.bash_profile`
  can take effect.
* You can verify that the installation worked by running:

```bash
conda --version
```

* Next, we will setup the the default channels (from where packages will
  be searched for and downloaded if no channel is specified).

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

# Practical exercise

## Get going with environments

Let's assume that you are just about to start a new exciting research project
called Project A.

* Let's make our first Conda environment:

```bash
conda create -n project_a -c bioconda fastqc
```

This will create an environment called project_a, containing FastQC from the
Bioconda channel. Conda will list the packages that will be installed and ask
for your confirmation.

* Once it is done, you can activate the environment:

```bash
conda activate project_a
```

By default, Conda will add information to your prompt telling you which
environment that is active.

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

* Did it work? Activate your project_a environment and run the `fastqc
  --version` command again. Does it work now?

Hopefully the FastQC software was not found in your base environment (unless
you had installed it previously), but worked once your environment was
activated.

* Now, let's add another package (sra-tools) to our environment using `conda
  install`. Make sure that project_a is the active environment first.

```bash
conda install -c bioconda sra-tools
```

* If we don't specify the package version, the latest available version will be
  installed. What version of sra-tools got installed?
* Run the following to see what versions are available:

```bash
conda search -c bioconda sra-tools
```

* Now try to install a different version of sra-tools, *e.g.*:

```bash
conda install -c bioconda sra-tools=2.7.0
```

Read the information that Conda displays in the terminal. It probably asks if
you want to downgrade the initial sra-tools installation to the one specified
here (2.7.0 in the example). You can only have one version of a given package
in a given environment.

Let's assume that you will have sequencing data in your Project A, and want to
use the latest bowtie2 software to align your reads.

* Find out what versions of bowtie2 are available in the bioconda channel using
  `conda search -c bioconda`.
* Now install the *latest* available version of bowtie2 in your project_a
  environment.

Let's further assume that you have an old project (called Project Old) where
you know you used bowtie2 v2.2.5. You just got back reviewer comments and they
want you to include some alignment statistics. Unfortunately, you haven't saved
that information so you will have to rerun the alignment. Now, it is essential
that you use the same version of bowtie that your results are based on,
otherwise the alignment statistics will be misleading. Using Conda environments
this becomes simple. You can just have a separate environment for your old
project where you have an old version of bowtie2 without interfering with your
new Project A where you want the latest version.

* Make a new environment for your old project:

```bash
conda create -n project_old -c bioconda bowtie2=2.2.5
```

* List your environments (do you remember the command?).
* Activate `project_old` and check the bowtie2 version (`bowtie2 --version`).
* Activate `project_a` again and check the bowtie2 version.
* Let's try to remove an installed package from the active environment:

```
conda remove sra-tools
```

* Run `conda deactivate` to exit your active environment.
* Now, let's remove an environment:

```bash
conda env remove -n project_old
```

After making a few different environments and installing a bunch of packages,
Conda can take up some disk space. You can remove unnecessary files with the
command:

```bash
conda clean -a
```

This will remove package tar-balls that are left from package installations,
unused packages (*i.e.* those not present in any environments), and cached
data.

!!! note "Quick recap"
    In this section we've learned:

    * How to use `conda install` for installing packages on the fly.
    * How to create and activate environments and how to change between them.
    * How to remove packages or environments and clean up.

## How to use in a reproducible project setting

We have up until now specified which Conda packages to install directly on the
command line using the `conda create` and `conda install` commands. For working
in projects this is not the recommended way. Instead, for increased control and
reproducibility, it is better to use a file (in [yaml
format](https://en.wikipedia.org/wiki/YAML)) specifying packages, versions and
channels needed to create the environment for a project.

* Throughout these tutorials we will use a case study where we analyze an
  RNA-seq experiment with the multiresistant bacteria MRSA (see
  [intro](tutorial_intro.md)). You will now start to make a Conda yml file for
  this MRSA project. The file will contain a list of the software and versions
  needed to execute the analysis code.
* In this conda tutorial, all code for the analysis is available in the script
  `code/run_qc.sh`. This code will:
    - download the raw fastq-files
    - run quality control on these using FastQC
* We will start by making a Conda yml-file that contains the required packages
  to perform these two steps.
* Later in the course, you will update the Conda yml-file with more packages,
  as the analysis workflow is expanded.
* So let's get going! Make a yaml file called `environment.yml` looking like
  this, and save it in the current directory (should be
  `reproducible_research_course/conda`):

```yaml
channels:
- conda-forge
- bioconda
dependencies:
- fastqc=0.11.9
- sra-tools=2.10.1
```

* Now, make a new Conda environment from the yaml file (note that here the
  command is `conda env create` as opposed to `conda create` that we used
  above):

```bash
conda env create -n project_mrsa -f environment.yml
```

!!! tip
    You can also specify exactly which channel a package should come from
    inside the environment file, using the `<channel>::<package>=<version>`
    syntax.

!!! tip
    Instead of the `-n` flag you can use the `-p` flag to set the full path to
    where the conda environment should be installed. In that way you can
    contain the conda environment inside the project directory, which does make
    sense from a reproducibility perspective, and makes it easier to keep track
    of what environment belongs to what project. If you don't specify `-p` the
    environment will be installed in the default `miniconda3/envs/` directory.

* Activate the environment!
* Now we can run the code for the MRSA project found in `code/run_qc.sh`,
  either by running `bash code/run_qc.sh` or by opening the `run_qc.sh` file
  and executing each line in the terminal one by one. Do this!

This should download the project fastq files and run FastQC on them (as
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

Projects can often be quite large and require lots of dependencies; it can feel daunting to try to capture all of that in a single Conda environment, especially when you consider potential incompatibilities that may arise. It can therefore be a good idea to start new projects with an environment file with each package you know that you will need to use, but without specifying exact version (except for those packages where you know you need a specific version). Conda will then try to get the latest compatible versions of all the specified software, making the start-up and installation part of new projects easier. You can then add the versions that were installed to your environment file afterwards, ensuring future reproducibility.

!!! note "Quick recap"
    In this section we've learned:

    * How to define our Conda environment using a yaml-file.
    * How to use `conda env create` to make a new environment from a yaml-file.

## Extra material

### Managing python versions

With conda it's possible to keep several different versions of python on your
computer at the same time, and switching between these versions is very easy. 
However, a single conda environment can only contain one version of python.

#### Your current python installation
The conda base environment has its own version of python installed.
When you open a terminal (after having installed conda on your system) this base
environment is activated by default (as evidenced by `(base)` prepended to your
prompt). You can check what python version is installed in this environment by 
running `python --version`. To see the exact path to the python executable type
`which python`.

In addition to this your computer may already have python installed in a 
separate (system-wide) location outside of the conda installation. To see if 
that is the case type `conda deactivate` until your prompt is not prepended
with a conda environment name. Then type `which python`. If a path was printed 
to the terminal (*e.g.* `/usr/bin/python`) that means some python version is
already installed in that location. Check what version it is by typing `python 
--version`.

Now activate the base conda environment again by typing `conda activate` (or
the equivalent `conda activate base`) then check the python installation path
and version using `which` and `python --version` as above. See the difference? 
When you activate a conda environment your `$PATH` variable is updated so that 
when you call `python` (or any other program) the system first searches the
directory of the currently active environment.

#### Different python versions

When you create a new conda environment you can choose to install a specific 
version of python in that environment as well. As an example, create an
environment containing python version `3.5` by running:

```bash
conda create -n py35 python=3.5
```

Here we name the environment `py35` but you can choose whatever name you want.

To activate the environment run:

```bash
conda activate py35
```

You now have a completely separate environment with its own python version. 

Let's say you instead want an environment with pythonversion `2.7` of installed.
You may for instance want to run scripts or packages that were written for
Python 2.x and are thus incompatible with Python 3.x. Simply create the new
conda environment with:

```bash
conda create -n py27 python=2.7
```

Activate this environment with:
```bash
conda activate py27
```

Now, switching between python versions is as easy as typing `conda activate 
py35` / `conda activate py27`.

!!! note "Default python version"
    If you create an environment where none of the packages require python, 
    **and** you don't explicitly install the `python` package then that new
    environment will use the python version installed in your base conda
    environment.

### Configuring conda

The behaviour of your conda installation can be changed using an optional
configuration file `.condarc`. On a fresh conda install no such file is
included but it's created in your home directory as `~/.condarc` the first time
you run `conda config`. 

You can edit the `.condarc` file either using a text editor or by way of the 
`conda config` command. To list all config parameters and their settings run:

```bash
conda config --show
```

Similar to conda environment files, the configuration file is in YAML syntax. 
This means that the config file is structured in the form of `key:value` pairs 
where the `key` is the name of the config parameter (*e.g.* `auto_update_conda`)
and the `value` is the parameter setting (*e.g.* `True`).

Adding the name of a config parameter to `conda config --show` will show only
that parameter, *e.g.* `conda config --show channels`.

You can change parameters with the `--set`, `--add`, `--append` and `--remove` 
flags to `conda config`. 

If you for example want to enable the 'Always yes' behaviour which makes conda 
automatically choose the `yes` option, such as when installing, you can run:

```bash
conda config --set always_yes True
``` 

To see details about a config parameter you can run `conda config --describe 
<parameter>`. Try running it on the `channels` parameter:

```bash
conda config --describe channels
```

In the beginning of this tutorial we added conda channels to the `.condarc`
file using `conda config --add channels`. To remove one of the channels from 
the configuration file you can run:

```bash
conda config --remove channels conda-forge
```

Check your `.condarc` file to see the change. To add the conda-forge channel
back to the top of the `channels` simply run:

```bash
conda config --add channels conda-forge
```

To completely remove a parameter and all its values run:

```bash
conda config --remove-key <parameter>
```

For a list of conda configuration parameters see the 
[Conda configuration](https://docs.conda.io/projects/conda/en/latest/configuration.html)
page.

## Tips and tricks

The following extra material contains some more advanced things you can do with
conda and the command line in general, which is not part of the main course
materials. All the essential skills of conda are covered by the previous
section: the material here should be considered tips and tricks from people who
use Conda as part of their daily work. You thus don't need to use these things
unless you want to, and you can even skip this part of the lesson if you like!

### Decorating your prompt

By default, Conda adds the name of the currently activated environment to the
end of your command line prompt. This is a good thing, as it makes it easier to
keep track of what environment and packages you have access to. The way this is
done in the default implementation becomes an issue when using absolute paths
for environments (specifying `conda env create -p <path/to/environment>`,
though, as the entire path will be added to the prompt. This can take up a lot
of unnecessary space, but can be solved in a number of ways.

The most straightforward way to solve this is to change the conda configuration
file, specifically the settings of the `env_prompt` configuration value which
determines how conda modifies your command line prompt. For more information
about this setting you can run `conda config --describe env_prompt` and to see
your current setting you can run `conda config --show env_prompt`.

By default env_prompt is set to '({default_env})' which modifies your prompt
with the active environment name if it was installed using the -n flag or if
the environment folder has a parent folder named envs/. Otherwise the full
environment path (*i.e.* the 'prefix') is displayed.

If you instead set env_prompt to '({name}) ' conda will modify your prompt with
the folder name of the active environment. You can change the setting by
running the following:

```bash
conda config --set env_prompt '({name}) '
```

If you wish to keep the ({default_env})' behaviour, or just don't want to
change your conda config, an alternative is to keep conda environment folders
within a parent folder called `envs/`. This will make conda only add the folder
name of the conda environment to your prompt when you activate it.

As an example, say you have a project called project_a with the project path
`~/myprojects/project_a`. You could then install the environment for project_a
into a folder `~/myprojects/project_a/envs/project_a_environment`. Activating
the environment by pointing conda to it (*e.g.* conda activate
`~/myprojects/project_a/envs/project_a_environment`) will only cause your promp
to be modified with project_a_environment.

### Bash aliases for conda

Some programmers like to have aliases (_i.e._ shortcuts) for common commands. 
Here are two aliases for conda that might prove useful for you.

```bash
# Activate an environment
alias coac='conda activate'

# Deactivate an environment
alias code='conda deactive'
```

Don't forget to add them to your `~/.bash_profile` if you want to use them!
