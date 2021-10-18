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
using `singularity build`. (Note that there is also the 
[singularity pull](https://sylabs.io/guides/3.4/user-guide/cli/singularity_pull.html)
command, which should be equivalent.)

Now, depending on the system you are running on and the version of Singularity, 
you may not have the option to build locally. However, Singularity has the 
option to build images remotely. To do this, you need to:

* Go to [https://cloud.sylabs.io/library](https://cloud.sylabs.io/library) and 
  create an account
* Log in and find "Access Tokens" in the menu and create a new token
* Copy the token
* In your terminal, run `singularity remote login` and hit ENTER. You should be
  asked to enter the token (API Key). Paste the copied token and hit ENTER. 
  You should get a **API Key Verified!** message.

> **Attention!** <br>
> In case you are not asked to enter the API Key, you can try to run 
> `singularity remote login SylabsCloud` instead.

We can now try to build the MRSA Singularity image using the `--remote` flag:

```bash
singularity build --remote mrsa_proj.sif docker://nbisweden/workshop-reproducible-research
```

This should result in a file called `mrsa_proj.sif`. 

In the Docker image we included the code needed for the workflow in the
`/course` directory of the image. These files are of course also available in
the Singularity image. However, a Singularity image is read-only (unless using
the sandbox feature), and this will be a problem if we try to run the workflow
within the `/course` directory, since the workflow will produce files and
Snakemake will create a `.snakemake` directory.  Instead, we need to provide
the files externally from our host system and simply use the Singularity image
as the environment to execute the workflow in (*i.e.* all the software).

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

## Containers at different platforms

A common problem with Singularity is that you can only create local builds if
you are working on a Linux system, as local builds for MacOS and Windows are
currently not supported. This means that you might favour using Docker instead
of Singularity, but what happens when you need to use a HPC cluster such as
Uppmax? Docker won't work there, as it requires root privileges, so Singularity
is the only solution. You can only run Singularity images there, however, not
*build* them...

So, how do you get a Singularity image for use on Uppmax if you can't build it
either locally or on Uppmax? You might think that using remote builds will
solve this, but for a lot of cases this won't help. Since most researchers will
want to work in private Git repositories they can't supply their Conda
`environment.yml` file to remote builds (which only works for public
repositories), which means that you’ll have to specify packages manually inside
the container instead.

There is, however, another solution: using Singularity inside Docker. By
creating a bare-bones, Linux-based Docker image with Singularity you can build
Singularity images locally on non-Linux operating systems. This can be either
from Singularity definition files or directly from already existing Docker
images. You can read more about this at the following [GitHub repository](https://github.com/kaczmarj/singularity-in-docker).
