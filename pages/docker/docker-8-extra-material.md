Docker is a large and complicated thing, but once you start using it regularly
you'll find that you start understand these complexities. There are lots of
different things you can do with images and containers in general, especially
when it comes to optimising build time or final image size. Here is some small
tips and tricks that you can be inspired from!

## A base image with Conda

We've used Conda throughout this Docker tutorial, and we did it by installing
Conda inside the image when we built it. Wouldn't it be nice if we didn't have
to do this particular step? After all, installing Conda is just busy-work,
compared to installing the actual environment that we want to use for the
analyses. Luckily, there are already Docker images out there that have Conda
(and [Mamba](conda-4-extra-material)) installed, such as the ones over at
`condaforge/mambaforge`! What follows is a Dockerfile that you could use instead
of the ones described above to install things using a Conda `environment.yml`
file, without having to install Conda in the Docker image when building it!

```no-highlight
FROM condaforge/mambaforge:4.10.1-0
LABEL description = "Image description"
MAINTAINER "Firstname Lastname" firstname.lastname@gmail.se

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /project

# Copy and install the Conda environment
COPY environment.yml ./
RUN conda config --set channel_priority strict \
    && mamba env update --name base --file environment.yml \
    && mamba clean --all --force-pkgs-dirs --yes

# Start Bash shell by default
CMD /bin/bash
```

## Singularity as an alternative container tool

Singularity is a container software alternative to Docker. It was originally
developed by researchers at Lawrence Berkeley National Laboratory with focus on
security, scientific software, and HPC clusters. One of the ways in which
Singularity is more suitable for HPC is that it very actively restricts
permissions so that you do not gain access to additional resources while inside
the container.

Here we give a brief introduction to Singularity and specifically how it can be
used on HPC clusters such as Uppmax.

If you want to read more, here are some additional resources:

* [Singularity docs](https://sylabs.io/guides/3.4/user-guide/index.html)
* [Uppmax Singularity user guide](
  https://www.uppmax.uu.se/support/user-guides/singularity-user-guide/)

### Converting Docker images to Singularity files

Singularity, unlike Docker, stores its images as a single file. A Singularity 
image file is self-contained (no shared layers) and can be moved around and 
shared like any other file.

While it is possible to define and build singularity images from scratch, in a
manner similar to what you've already learned for Docker, that is not something
we will cover here (but feel free to read more about this in _e.g._ the 
[Singularity docs](https://sylabs.io/guides/3.4/user-guide/quick_start.html#singularity-definition-files)).

Instead, we will take advantage of the fact that Singularity can convert Docker 
images to the Singularity Image Format (SIF). This is great if there's a Docker 
image that you want to use on an HPC cluster such as Uppmax where you cannot use 
Docker. 

Let's try to convert the Docker image for this course directly from DockerHub 
using Singularity:

```bash
singularity pull mrsa_proj.sif docker://nbisweden/workshop-reproducible-research
```

This should result in a file called `mrsa_proj.sif`. 

In the Docker image we included the code needed for the workflow in the
`/course` directory of the image. These files are of course also available
in the Singularity image. However, a Singularity image is read-only (unless
using the sandbox feature), and this will be a problem if we try to run the
workflow within the `/course` directory, since the workflow will produce
files and Snakemake will create a `.snakemake` directory. 

vide the files externally from our host system and
simply use the Singularity image as the environment to execute the workflow
in (i.e. all the software). 

In your current working directory (`workshop-reproducible-research/tutorials/docker/`) 
the vital MRSA project files are already available (`Snakefile`, `config.yml`, 
`code/header.tex` and `code/supplementary_material.Rmd`). 

Since Singularity bind mounts the current working directory we can simply
execute the workflow and generate the output files using:

```bash
singularity run --vm-ram 2048 mrsa_proj.sif
```

This executes the default run command, which is 
`snakemake -rp --configfile config.yml` (as defined in the original 
`Dockerfile`). 

> **Note** <br>
> Note here that we have increased the allocated RAM to 2048 MiB (`--vm-ram 2048`), 
> needed to fully run through the workflow. In case the command fails, 
> you can try to increase the RAM to e.g. 4096 MiB, or you can try to run the
> command without the  `--vm-ram` parameter.

The previous step in this tutorial included running the `run_qc.sh` script, 
so that part of the workflow has already been run and Snakemake will continue 
from that automatically without redoing anything. Once completed you should 
see a bunch of directories and files generated in your current working 
directory, including the `results/` directory containing the final HTML report.