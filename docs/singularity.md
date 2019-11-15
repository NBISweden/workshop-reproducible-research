# Introduction to Singularity

## What is Singularity?

## Tell me more

* [Singularity docs](https://sylabs.io/guides/3.4/user-guide/index.html)
* [Uppmax Singularity user guide](https://www.uppmax.uu.se/support/user-guides/singularity-user-guide/)

# Setup

This tutorial depends on files from the course Bitbucket repo. Take a look at the [intro](tutorial_intro.md) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/singularity`.

## Install Singularity

### macOS

Download the Singularity Desktop DMG file from [here](https://sylabs.io/singularity-desktop-macos/) and follow the instructions. Note that this is a beta version and not all features are available yet.

### Linux

Follow the instructions [here](https://sylabs.io/guides/3.4/user-guide/installation.html#distribution-packages-of-singularity).

### Windows

Installing on Windows may be tricky. This requires running Singularity through a Vagrant Box. See [instructions here](https://sylabs.io/guides/3.4/user-guide/installation.html#install-on-windows-or-mac).

# Practical exercise

## The very basics

In the [Docker tutorial](docker.md) we started by downloading an Ubuntu image. Let's see how the same thing is achieved with Singularity:

```bash
singularity pull library://ubuntu
```

This pulls an ubuntu image from the [Singularity library](https://cloud.sylabs.io/library) (somewhat equivalent to Dockerhub). The first thing you might have noticed is that this command produces a file `ubuntu_latest.sif` in the current working directory. Singularity, unlike Docker, stores its images as a single file. Docker on the other hand uses layers, which can be shared between multiple images, and thus stores downloaded images centrally (remember the `docker image ls` command?). A Singularity image file is self-contained (no shared layers) and can be moved around and shared like any other file.

To run a command in a Singularity container (equivalent of e.g. `docker run ubuntu uname -a`) we can execute:

```bash
$ singularity exec ubuntu_latest.sif uname -a
Linux (none) 4.19.10 #1 SMP Mon Apr 8 00:07:40 CDT 2019 x86_64 x86_64 x86_64 GNU/Linux
[    4.994162] reboot: Power down
```

Now, try to also run the following commands in the ubuntu container in the same manner as above:

* `whoami`
* `ls -lh`

Notice anything unexpected or different from what you learnt from the Docker tutorial?

Unlike Docker, Singularity attempts to map parts of your local file system to the image. By default Singularity bind mounts `$HOME`, `/tmp`, and `$PWD` (the current working directory) into your container. Also, inside a Singularity container, you are the same user as you are on the host system.

We can also start an interactive shell (equivalent of e.g. `docker run -it ubuntu`):

```bash
singularity shell ubuntu_latest.sif
```

While running a shell in the container, try executing `pwd` (showing the full path to your current working directory). See that it appears to be your local working directory? Try `ls /` to see the files and directory in the root. Exit the container (`exit`) and run `ls /` to see how it looks on your local system. Notice the difference?

!!! note "Quick recap"
    In this section we covered:

    * how to download a Singularity image using `singularity pull`
    * how to run a command in a Singularity container using `singularity exec`
    * how to start an interactive terminal in a Singularity container using `singularity shell`


## Bind mounts

In the previous section we saw how Singularity differs from Docker in terms of images being stored in stand-alone files and much of the host filesystem being mounted in the container. We will now explore this further.

Similarly to the `-v` flag for Docker we can use `-B` or `--bind` to bind mount directories into the container. For example, to mount a directory into the `/mnt/` directory in the container one would do:

```bash
singularity shell -B /path/to/dir:/mnt ubuntu_latest.sif
```

You can try this for example by mounting the `conda/` tutorial directory to `/mnt`:

```bash
singularity shell -B ../conda:/mnt ubuntu_latest.sif
```

In the container, to see that the bind mounting worked, run e.g.:

```bash
ls /mnt/code
```

Now, this was not really necessary since `conda/` would have been available to us anyway since it most likely has you home directory as a parent, but it illustrates the capabilities to get files from the host system into the container when needed. Note also that if you run Singularity on say an HPC cluster, the system admins may have enabled additional default directories that are bind mounted automatically.

!!! note "Quick recap"
    In this section we covered:

    * how to bind mount specific directories using `-B`

## Pulling Docker images

Singularity has the ability to convert Docker images to the Singularity Image Format (SIF). We can try this out by running:

```bash
singularity pull docker://godlovedc/lolcow
```

This command generates a .sif file where the individual layers of the specified Docker image have been combined and converted to Singularity's native format. We can now use `run`, `exec`, and `shell` commands on this image file. Try it:

```bash
singularity run lolcow_latest.sif
singularity exec lolcow_latest.sif fortune
singularity shell lolcow_latest.sif
```



!!! note "Quick recap"
    In this section we covered:

    * how to use `singularity pull` to download and run Docker images as Singularity containers

## Building a Singularity image from scratch

As we have seen, it is possible to convert Docker images to the Singularity format when needed and run them using Singularity. In terms of making a research project reproducible using containers, it may be enough to e.g. define a Dockerfile (recipe for a Docker image) as well as supply a Docker image for others to download and use, either directly through Docker, or by Singularity. Even better, from a reproducibility aspect, would be to also generate the Singularity image from the Docker image and provide that for potential future users (since the image is a static file, whereas running `singularity pull` or `singularity build` would rebuild the image at the time of issuing the command).

A third option, would be to define a Singularity recipe, either on its own or in addition to the Dockerfile. The equivalent of a Dockerfile for Singularity is called a Singularity Definition file ("def file").

The def file consists of two parts:

* a header that defines the core OS and related features
* optional sections, each starting with a `%`-sign, that add content or execute commands during the build process

As an example, let's look at the def file for the lol_cow image that we played with above (above we pulled lol_cow from Dockerhub but it also exists in the Singularity library and can be pulled by `singularity pull library://godlovedc/funny/lolcow`).

The file looks like this:

```
BootStrap: docker
From: ubuntu:16.04

%post
    apt-get -y update
    apt-get -y install fortune cowsay lolcat

%environment
    export LC_ALL=C
    export PATH=/usr/games:$PATH

%runscript
    fortune | cowsay | lolcat
```

The first part of the header sets the bootstrap agent. In the lol_cow example DockerHub is used. Alternatively one could set it to *library* to use the Singularity Library. There are also other bootstrap agents available (see [this link](https://sylabs.io/guides/3.3/user-guide/definition_files.html#preferred-bootstrap-agents) for details). Next, the base image that the new image starts from is defined, in this case the Docker image `ubuntu:16.04`.

In the lol_cow def file three sections are used (`%post`, `%environment`, and `runscript`).

* `%post` is similar to the `RUN` instruction in Dockerfiles. Here is where you include code to download files from the internet, install software, create directories etc.
* `%environment` is similar to the `ENV` instruction in Dockerfiles. It is used to set environmental variables that will be available when running the container. The variables set in this section will not however be available during the build and should in the cases they are needed then also be set in the `%post` section.
* `%runscript` is similar to the `CMD` instruction in Dockerfiles and contains the default command to be executed when running the container.

### Singularity def file for the MRSA project

Let's use the MRSA case study project to define our own Singularity def file! We will not make an image for the whole workflow but rather focus on the `run_qc.sh` script that we used in the end of the [conda tutorial](conda.md). This script is included in the `code` directory of your current working directory (`singularity`) and, when executed, downloads a few fastq-files and runs FastQC on them. To run the script we need the software SRA-Tools and FastQC.

* Make a new file called `run_qc.def` and add the following lines:

```
Bootstrap: library
From: ubuntu:18.04
```

Here we'll use the Singularity Library as bootstrap agent, instead of Dockerhub as in the lol_cow example above. The base Singularity image will be `ubuntu:18.04`.

* Next, add the `%environment` section:

```
%environment
    export LC_ALL=C
    export PATH=/usr/miniconda3/bin:$PATH
```

This sets the default locale as well as includes the PATH to conda (which we will soon install).

%files
%test
%labels

* Now add the `%post` section:

```
%post
    apt-get update
    apt-get install -y --no-install-recommends bzip2 ca-certificates curl
    apt-get clean

    # Install conda:
    curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O
    bash Miniconda3-4.6.14-Linux-x86_64.sh -bf -p /usr/miniconda3/
    rm Miniconda3-4.6.14-Linux-x86_64.sh

    # Configure conda:
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # Install requirements:
    conda install -c bioconda fastqc=0.11.6 sra-tools=2.8
    conda clean --all
```

You should recognize parts of this from the Docker tutorial. Basically, we install some required basic tools like bzip2, cca-certificates and curl, then install and configure conda and finally install the required tools for the `run_qc.sh` script.

* Next, add a `%test` section:

```
%test
    fastqc --version
    fastq-dump --version
```

The test section runs at the end of the build process and you can include any code here to verify that your image works as intended. Here we just make sure that the `fastqc` and `fastq-dump` are available.

Finally, add the `%runscript`:

```
%runscript
    bash code/run_qc.sh
```

We should now be ready to build our image from this def file. Do that by running:

```bash
singularity build --remote run_qc.sif run_qc.def
```

Did it work? Can you figure out why it failed? Tip: it has to do with the PATH.

??? note "Click to see the solution"
    You need to add conda to the PATH. `%environment` makes it available at runtime but not during build.

    Update the `%post` section as follows:

    ```bash
    # Install conda:
    curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O
    bash Miniconda3-4.6.14-Linux-x86_64.sh -bf -p /usr/miniconda3/
    rm Miniconda3-4.6.14-Linux-x86_64.sh
    export PATH=/usr/miniconda3/bin:$PATH ## <- add this line
    ```

    You also need to update the `%test` section:

    ```bash
    export PATH=/usr/miniconda3/bin:$PATH ## <- add this line
    fastqc --version
    fastq-dump --version
    ```

The build should now hopefully work and produce a Singularity image called `run_qc.sif`. To run the image, thus executing `code/run_qc.sh` using the tools in the container, do:

```bash
singularity run run_qc.sif
```

!!! tip
    For building testing convenience one can use sandbox, shell and install, then add it to he def file for final build


## Converting the MRSA workflow Docker image to Singularity

Instead of `pull` one can use `build` to build Singularity images with access to more options, e.g. like building a writable sandbox image. We will use `singularity build` to convert the Docker image of the MRSA project, that we use as a case study in this course, to a Singularity image. Now depending on the system you are running on and the version of Singularity you may not have the option to build locally. However, Singularity has the option to build images remotely. To do this, you need to:

* Go to [https://cloud.sylabs.io/library](https://cloud.sylabs.io/library) and create an account
* Log in and find "Access Tokens" in the menu and create a new token
* Copy the token and paste it in the file `~/.singularity/sylabs-token`

We can now try to build the MRSA Singularity image using the `--remote` flag:

```bash
singularity build --remote mrsa_proj.sif docker://scilifelablts/reproducible_research_course
```

This should result in a file called `mrsa_proj.sif`. In the Docker image we included the code needed for the workflow in the `/course` directory of the image. These files are of course also available in the Singularity image. However, a Singularity image is read-only (unless using the sandbox feature), and this will be a problem if we try to run the wworkflow within the `/course` directory, since the workflow will produce files and Snakemake will create a `.snakemake` directory. Instead, we need to provide the files externally from our host system and simply use the sSingularity image as the environment to execute the workflow in (i.e. all the software). In your current working directory (`singularity/`) the vital MRSA project files are already available (`Snakefile`, `config.yml`, `code/header.tex` and `code/supplementary_material.Rmd`). And since Singularity bind mounts the current working directory we can simply execute the workflow and generate the output files using:

```bash
singularity run --vm-ram 2048 mrsa_proj.sif
```

This executes the default run command, which is `snakemake -rp --configfile config.yml` (as defined in the original `Dockerfile`). Note here that we have also increased the allocated RAM to 2048 MiB (`--vm-ram 2048`), needed to fully run through the workflow. Once completed you should see a bunch of directories and files generated in your current working directory, including the `results/` directory containing the final PDF report.

### Other
inspect
instances
more options to shell run exec
