# Introduction to Singularity

## What is Singularity?

## Tell me more

* [Singularity docs](https://sylabs.io/guides/3.4/user-guide/index.html)

# Setup

This tutorial depends on files from the course Bitbucket repo. Take a look at the [intro](tutorial_intro.md) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/singularity`.

## Install Singularity

### macOS

Download the Singularity Desktop DMG file from [here](https://sylabs.io/singularity-desktop-macos/) and follow the instructions.

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

This pulls an ubuntu image from the [Singularity library](https://cloud.sylabs.io/library). The first thing you might have noticed is that this command produces a file `ubuntu_latest.sif` in the current working directory. Singularity, unlike Docker, stores its images as a single file. Docker on the other hand uses layers, which can be shared between multiple images, and thus stores downloaded images centrally (remember the `docker image ls` command?). A Singularity image file is self-contained (no shared layers) and can be moved around and shared like any other file.

To run a command in a Singularity container (equivalent of e.g. `docker run ubuntu uname -a`) we can execute:

```bash
$ singularity exec ubuntu_latest.sif uname -a
Linux (none) 4.19.10 #1 SMP Mon Apr 8 00:07:40 CDT 2019 x86_64 x86_64 x86_64 GNU/Linux
[    4.994162] reboot: Power down
```

Now, try to also run the commands `whoami` and `ls -lh` in the ubuntu container in the same manner as above. Notice anything unexpected or different from what you learnt from the Docker tutorial?

Unlike Docker, Singularity attempts to map parts of your local file system to the image. User. Default mapped dirs. Env vars.

We can also start an interactive shell (equivalent of e.g. `docker run -it ubuntu`):

```bash
singularity shell ubuntu_latest.sif
```

While running a shell in the container, try executing `pwd` (showing the full path to your current working directory). See that it appears to be your local working directory? Try `ls /` to see the files and directory in the root. Exit the container and run `ls /` to see how it looks on your local system. Notice the difference?
