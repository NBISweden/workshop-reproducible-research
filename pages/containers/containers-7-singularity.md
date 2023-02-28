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

* [Singularity docs](https://sylabs.io/guides/master/user-guide/)
* [Uppmax Singularity user guide](
  https://www.uppmax.uu.se/support/user-guides/singularity-user-guide/)

!!! Info "Singularity and Apptainer"
    Singularity has very recently been renamed to *Apptainer*, but we have opted
    to stick with the original name in the material for now, while the change is
    still being adopted by the community and various documentation online.

### Converting Docker images to Singularity files

Singularity, unlike Docker, stores images as single files. A Singularity
image file is self-contained (no shared layers) and can be moved around and
shared like any other file.

While it is possible to define and build Singularity images from scratch, in a
manner similar to what you've already learned for Docker, this is not something
we will cover here (but feel free to read more about this in _e.g._ the
[Singularity docs](https://sylabs.io/guides/master/user-guide/definition_files.html)).

Instead, we will take advantage of the fact that Singularity can convert Docker
images to the Singularity Image Format (SIF). This is great if there's a Docker
image that you want to use on an HPC cluster such as Uppmax where you cannot use
Docker.

!!! Tip
    If you are running singularity through Vagrant VirtualBox you may have to
    set the temporary directory that Singularity uses during pull/build commands
    to something with more disk space. First run `mkdir ~/tmp` to create a tmp
    directory inside the home folder of the VirtualBox, then
    `export SINGULARITY_TMPDIR="~/tmp"`.

Let's try to convert the Docker image for this course directly from DockerHub
using `singularity pull`:

```bash
singularity pull mrsa_proj.sif docker://nbisweden/workshop-reproducible-research
```

This should result in a file called `mrsa_proj.sif`.

### Running a singularity image

In the Docker image we included the code needed for the workflow in the
`/course` directory of the image. These files are of course also available in
the Singularity image. However, a Singularity image is read-only (unless using
the [sandbox](https://sylabs.io/guides/master/user-guide/build_a_container.html#creating-writable-sandbox-directories)
feature). This will be a problem if we try to run the workflow
within the `/course` directory, since the workflow will produce files and
Snakemake will create a `.snakemake` directory.  Instead, we need to provide
the files externally from our host system and simply use the Singularity image
as the environment to execute the workflow in (*i.e.* all the software and
dependencies).

In your current working directory (`workshop-reproducible-research/tutorials/containers/`)
the vital MRSA project files are already available (`Snakefile`, `config.yml`,
`code/header.tex` and `code/supplementary_material.Rmd`).

Since Singularity bind mounts the current working directory we can simply
execute the workflow and generate the output files using:

```bash
singularity run mrsa_proj.sif
```

This executes the default run command, which is
`snakemake -rp -c 1 --configfile config.yml` (as defined in the original
`Dockerfile`).

The previous step in this tutorial included running the `run_qc.sh` script,
so that part of the workflow has already been run and Snakemake will continue
from that automatically without redoing anything. Once completed you should
see a bunch of directories and files generated in your current working
directory, including the `results/` directory containing the final HTML report.

## Containers at different platforms

A common problem with Singularity is that you can only create local builds if
you are working on a Linux system, as local builds for MacOS and Windows are
currently not supported. This means that you might favour using Docker instead
of Singularity, but what happens when you need to use a HPC cluster such as
Uppmax? Docker won't work there, as it requires root privileges, so Singularity
is the only solution. You can only run Singularity images there, however, not
*build* them...

So, how do you get a Singularity image for use on Uppmax if you can't build it
either locally or on Uppmax? While it's possible to do remote builds (via the
`--remote` flag), in our experience this functionality is not stable and for a
lot of cases it won't help. Since most researchers will want to work in private
Git repositories they can't supply their Conda `environment.yml` file to remote
builds (which only works for public repositories), which means that you’ll have
to specify packages manually inside the container instead.

There is, however, another solution: using Singularity inside Docker. By
creating a bare-bones, Linux-based Docker image with Singularity you can build
Singularity images locally on non-Linux operating systems. This can be either
done from Singularity definition files or directly from already existing Docker
images. You can read more about this at the following [GitHub repository](https://github.com/kaczmarj/singularity-in-docker).


!!! Success "Quick recap"
    In this section we've learned:

    - How to convert Docker images to Singularity images.
    - How to use `singularity run` for starting a container from an image.
    - How to build a Singularity image using Singularity inside Docker.
