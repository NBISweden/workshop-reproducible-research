# Introduction to Singularity

## What is Singularity?

## Tell me more

* [Singularity docs](https://sylabs.io/guides/3.4/user-guide/index.html)

# Setup

This tutorial depends on files from the course Bitbucket repo. Take a look at the [intro](tutorial_intro.md) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/singularity`.

## Install Singularity

### macOS

Download the Singularity Desktop DMG file from [here](https://sylabs.io/singularity-desktop-macos/) and follow the instructions. Note that this is a beta version and not all features are available yet.

### Linux

Follow the instructions [here](https://sylabs.io/guides/3.4/user-guide/installation.html#distribution-packages-of-singularity).

### Windows

Installing on Windows may be tricky. This requires running Singularity through a Vagrant Box. See [instructions here](https://sylabs.io/guides/3.4/user-guide/installation.html#install-on-windows-or-mac).

## Practical exercise

### The very basics

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


### Bind mounts

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

### Pulling Docker images

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

Instead of `pull` one can use `build` to build Singularity images with access to more options, e.g. like building a writable sandbox image. We will use `singularity build` to convert the Docker image of the MRSA project, that we use as a case study in this course, to a Singularity image. Now depending on the system you are running on and the version of Singularity you may not have the option to build locally. However, Singularity has the option to build images remotely. To do this, you need to:

* Go to https://cloud.sylabs.io/library and create an account
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

This executes the default run command, which is `snakemake -rp --configfile config.yml`. Note here that we have also increased the allocated RAM to 2048 MiB (`--vm-ram 2048`), needed to fully run through the workflow. Once completed you should see a bunch of directories and files generated in your current working directory, including the `results/` directory containing the final PDF report.

### Building a Singularity image

definition file
