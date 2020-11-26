# Introduction

Welcome to the tutorials! Here we will learn how to make a computational
research project reproducible using several different tools, described in the
figure below:

![](images/tutorials_overview.png)

The figure gives an overview of the six available tutorials, a very brief
description of their main purpose, and the suggested order to do them. However,
each tutorial is made so that it can be completed independently of the other
tutorials. It is therefore perfectly possible to choose a different order, or
a subset of tutorials that suits your interests. Under the main figure there is
a list of a few suggested alternative tutorial orders; you will find the
tutorials in the menu to the left.

Before going into the tutorials themselves, we first describe the case study
from which the example data comes from, followed by the setup needed to install
the tools themselves. These will create quite a lot of files on your computer,
some of which will actually take up a bit of storage space too. In order to
remove any traces of these after completing the tutorials, please refer to the
[Take down section](take_down.md).

## The case study

We will be running a small bioinformatics project as a case study, and use that
to exemplify the different steps of setting up a reproducible research project.
To give you some context, the study background and analysis steps are briefly
described below.

### Background

The data is taken from [Osmundson, Dewell, and Darst (2013)](
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076572),
who have studied methicillin-resistant *Staphylococcus aureus* (MRSA). MRSA is
resistant to broad spectrum beta-lactam antibiotics and lead to
difficult-to-treat infections in humans. Lytic bacteriophages have been
suggested as potential therapeutic agents, or as the source of novel antibiotic
proteins or peptides. One such protein, gp67, was identified as
a transcription-inhibiting transcription factor with an antimicrobial effect.
To identify *S. aureus* genes repressed by gp67, the authors expressed gp67 in
*S. aureus* cells. RNA-seq was then performed on three S. aureus strains:

* RN4220 with pRMC2 with gp67
* RN4220 with empty pRMC2
* NCTC8325-4

### Analysis

The graph below shows the different steps of the analysis that are included in
this project:

![](images/rulegraph_mrsa_intro.svg)

The input files are:

* RNA-seq raw data (FASTQ files) for the three strains
* *S. aureus* genome sequence (a FASTA file)
* *S. aureus* genome annotation (a GFF file)

The different steps of the workflow and what they do are as follows:

* `get_genome_fasta` - Download the genome file.
* `index_genome` - Index the genome using the *Bowtie2* software (required for
  the alignment step)
* `get_SRA_by_accession` - Download the RNA-seq raw data for the three strains
  from the *Sequence Read Archive* (SRA).
* `fastqc` - Run quality control on each of the RNA-seq FASTQ files using the
  *FastQC* software.
* `multiqc` - Summarize the quality controls.
* `align_to_genome` - Align the RNA-seq data from the three strains to the
  indexed genome using the *Bowtie2* software.
* `sort_bam` - Sort the alignment files by genome coordinate using the
  *Samtools* software.
* `get_genome_gff3` - Download the genome annotation file.
* `generate_count_table` - Calculate gene expression by counting aligned reads
  per gene using the *HTSeq-count* software.
* `generate_rulegraph` - Generate the workflow overview figure shown above.
* `make_supplementary` - Produce the supplementary materials section using data
  from the quality controls, gene counts and the workflow figure.

## Setup for Mac / Linux users

Clone the GitHub repository containing all files you will need for completing
the tutorials. First, `cd` into a directory on your computer (or create one)
where it makes sense to download the course directory.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

!!! attention
    Check your git version with `git --version`. If you have a very old version
    of git you might want to update to a later version.

!!! tip
    If you want to revisit the material from an older instance of this course,
    you can do that using `git checkout tags/<tag-name>`, e.g. `git checkout
    tags/course_1905`. To list all available tags, use `git tag`. Run this
    command after you have `cd` into `workshop-reproducible-research` as
    described above. If you do that, you probably also want to view the same
    older version of this website. Locate the version box in the bottom right
    corner of the browser and select the corresponding version.

## Setup for Windows users

There are several different ways to run the course material on a Windows
computer. Neither is perhaps optimal, and the material itself has not been
adapted specifically for Windows. Nevertheless, in principle everything
*should* be possible to run. A few ways you could setup:

- Use the Linux Bash Shell on Windows 10 (see below) _**Recommended for the
  course**_
- Run the course in a Docker container (see below)
- Run as Linux through a virtual machine (and see the Linux setup above)
- Use the Windows 10 PowerShell, install git and clone the course repository
  (see the Mac/Linux setup above)

### Running in the Linux Bash Shell on Windows 10

This will give you access to a full command-line bash shell based on Linux on your
Windows 10 PC. For the difference between the Linux Bash Shell and the PowerShell on Windows
10, see *e.g.* [this article](
https://searchitoperations.techtarget.com/tip/On-Windows-PowerShell-vs-Bash-comparison-gets-interesting).

Install Bash on Windows 10, following the instructions at *e.g.* one of these
resources:

- [Installing the Windows Subsystem and the Linux Bash](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
- [Installing and using Linux Bash on Windows](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
- [Installing Linux Bash on Windows](https://itsfoss.com/install-bash-on-windows/)

Open a bash shell Linux terminal and clone the GitHub repository 
containing all files you will need for completing the tutorials as follows. 
First, `cd` into a directory on your computer (or create one) where it makes 
sense to download the course directory.

!!! tip
    You can find the directory where the Linux distribution is storing all its files by
    typing `explorer.exe .`. This will launch the Windows File Explorer showing the 
    current Linux directory.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

!!! tip
    If you want to revisit the material from an older instance of this course,
    you can do that using `git checkout tags/<tag-name>`, *e.g.* `git checkout
    tags/course_1905`. To list all available tags, use `git tag`. Run this
    command after you have `cd` into `workshop-reproducible-research` as
    described above. If you do that, you probably also want to view the same
    older version of this website. Locate the version box in the bottom right
    corner of the browser and select the corresponding version.

### Using Docker to run the course

Alternatively, you can use Docker to run the course in a Docker container.
First, open the Windows 10 PowerShell and `cd` into a directory on your computer 
(or create one) where it makes sense to download the course directory. 
Install Docker by following the instructions in the [Docker tutorial](docker.md#windows). 
Then run:

```bash
cd c:/my_dir
docker run -it -p 8888:8888 -v /c/my_dir:/course/ \
    nbisweden/workshop-reproducible-research:slim
```

!!! attention
    Note that we use `/c/my_dir` and not `c:/my_dir` as we normally do on
    Windows. This is required for Docker to parse the command correctly.

This will start an isolated container running Linux, where the directory
`c:/my_dir` is mounted (*i.e.* you can access the files in this Windows
directory within the Linux container, and files edited or created within the
Linux container will appear in this Windows directory). Note that the idea is
that you should edit files in the mounted `c:/my_dir` using an editor in your
normal OS, say Notepad in Windows. The terminal in the container is for running
stuff, not editing.

You should now be at a terminal in the Docker container. Now clone the GitHub
repository containing all the files you will need for the tutorials.

```bash
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

#### Setting up persistent conda environments
Because the root file system of each container is an isolated instance conda
environments you create during this course will be lost if you exit the 
container or if it is killed for some reason. This means that you will have to 
recreate each environment every time you run a new container for
the course. To avoid this, you can make sure conda uses a subfolder `envs/` 
inside the `/course` directory for storing environments. If you run containers
for the course with some folder on your local machine mounted inside `/course`
that will cause the conda environments to be available on your local machine
even though they are created inside the container. Should a container be stopped
for some reason you can simply run a new one and activate conda environments
under `/course/envs`, saving you the trouble of recreating them.

What you have to do is to, after you've started a container with the 
`docker run` command above, first create the `envs/` directory:

```bash
mkdir /course/envs
``` 

Then set the environment variable `CONDA_ENVS_PATH` to `/course/envs`:

```bash
export CONDA_ENVS_PATH="/course/envs"
```

Now when you use Conda to create environments inside the Docker container for
the course, the environments will be placed in `/course/envs` which is also 
present on your local system because you mounted the folder inside the 
container. **Note that you will have to set the `CONDA_ENVS_PATH` each time you
start a new container for this to work**.

!!! attention
    This usage of persistent conda environments should be considered an edge 
    case of how you use Docker containers. We do this only to make it easier to
    run the course through Docker, but in general we do not advocate creating
    conda environments separate from the actual Docker container.  

Don't worry if you feel that this Docker stuff is a little confusing, it will
become clearer in the [Docker tutorial](docker.md). However, the priority right
now is just to get it running so that you can start working.
