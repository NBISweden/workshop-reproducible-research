## Apptainer as an alternative container tool

Apptainer is a container software alternative to Docker. It was originally
developed as _Singularity_ by researchers at Lawrence Berkeley National
Laboratory (read more about this below) with focus on security, scientific
software, and HPC clusters. One of the ways in which Apptainer is more suitable
for HPC is that it very actively restricts permissions so that you do not gain
access to additional resources while inside the container. Apptainer also,
unlike Docker, stores images as single files using the _Singularity Image
Format_ (SIF). A SIF file is self-contained and can be moved around and shared
like any other file, which also makes it easy to work with on an HPC cluster.

> **Singularity and Apptainer** <br>
> The open source Singularity project was renamed to _Apptainer_ in 2021.
> Confusingly, the company _Sylabs_ still keeps their commercial branch of
> the project under the Singularity name, and offer a free 'Community
> Edition' version. The name change was done in order to clarify the
> distinction between the open source project and the various commercial
> versions. At the moment there is virtually no difference to you as a user
> whether you use Singularity or Apptainer, but eventually it's very likely that
> the two will diverge.

While it is possible to define and build Apptainer images from scratch, in a
manner similar to what you've already learned for Docker, this is not something
we will cover here (but feel free to read more about this in _e.g._ the
[Apptainer docs](https://apptainer.org/docs/user/main/index.html).

The reasons for not covering Apptainer more in-depth are varied, but it
basically boils down to it being more or less Linux-only, unless you use Virtual
Machines (VMs). Even with this you'll run into issues of incompatibility of
various kinds, and these issues are further compounded if you're on one of the
new ARM64-Macs. You also need `root` (admin) access in order to actually _build_
Apptainer images regardless of platform, meaning that you can't build them on
_e.g._ Uppmax, even though Apptainer is already installed there. You can,
however, use the `--remote` flag, which runs the build on Apptainer's own
servers. This doesn't work in practice a lot of the time, though, since most
scientist will work in private Git repositories so that their research and code
is not available to anybody, and the `--remote` flag requires that _e.g._ the
`environment.yml` file is publicly available.

There are very good reasons to use Apptainer, however, the major one being
that you aren't allowed to use Docker on most HPC systems! One of the nicer
features of Apptainer is that it can convert Docker images directly for use
within Apptainer, which is highly useful for the cases when you already built
your Docker image or if you're using a remotely available image stored on _e.g._
DockerHub. For a lot of scientific work based in R and/or Python, however, it is
most often the case that you build your own images, since you have a complex
dependency tree of software packages not readily available in existing images.
So, we now have another problem for building our own images:

1. Only Apptainer is allowed on HPC systems, but you can't build images there
   due to not having `root` access.
2. You can build Apptainer images locally and transfer them to HPCs, but this
   is problematic unless you're running Linux natively.

Seems like a "catch 22"-problem, right? There are certainly workarounds (some of
which we have already mentioned) but most are roundabout or difficult to get
working for all use-cases. Funnily enough, there's a simple solution: run
Apptainer locally from inside a Docker container! Conceptually very meta, yes,
but works very well in practice. What we are basically advocating for is that
you stick with Docker for most of your container-based work, but convert your
Docker images using Apptainer-in-Docker whenever you need to work on an HPC.
This is of course not applicable to Linux users or those of you who are fine
with working through using VMs and managing any issues that arise from doing
that.

> **Summary** <br>
> Apptainer is a great piece of software that is easiest to use if you're
> working on a Linux environment. Docker is, however, easier to use from a
> cross-platform standpoint and covers all use-cases except running on HPCs.
> Running on HPCs can be done by converting existing Docker images at runtime,
> while building images for use on HPCs can be done using local Docker images
> and Apptainer-in-Docker.

## Apptainer-in-Docker

By creating a bare-bones, Linux-based Docker image with Apptainer you can
build Apptainer images locally on non-Linux operating systems. There is
already a good image setup for just this, and it is defined in this [GitHub
repository](https://github.com/kaczmarj/apptainer-in-docker). Looking at the
instructions there we can see that we need to do the following:

```bash
docker run \
    --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v $(pwd):/work \
    kaczmarj/apptainer \
    build <IMAGE>.sif docker-daemon://<IMAGE>:<TAG>
```

You already know about `docker run`, the `--rm` flag and bind mounts using `-v`.
The `/var/run/docker.sock` part is the Unix socket that the Docker daemon
listens to by default, meaning that it is needed for us to be able to
specify the location of the Docker container we want to convert to a SIF
file. The `kaczmarj/apptainer` part after the bind mounts is the image
location hosted at [DockerHub](https://hub.docker.com/r/kaczmarj/apptainer),
while the last line is the Apptainer command that actually does the conversion.
All we need to do is to replace the `<IMAGE>` part with the Docker image we want
to convert, _e.g._ `my_docker_image`.

- Replace `<IMAGE>` and `<TAG>` with one of your locally available Docker images
  and one of its tags and run the command - remember that you can use `docker
image ls` to check what images you have available.

In the end you'll have a SIF file (_e.g._ `my_docker_image.sif`) that you can
transfer to an HPC such as Uppmax and run whatever analyses you need. If you
want to be able to do this without having to remember all the code you can check
out the [this script](https://github.com/fasterius/dotfiles/blob/main/scripts/apptainer-in-docker.sh).

## Running Apptainer

The following exercises assume that you have a login to the Uppmax HPC cluster
in Uppsala, but will also work for any other system that has Apptainer
installed - like if you managed to install Apptainer on your local system or
have access to some other HPC cluster. Let's try to convert the Docker image for
this course directly from DockerHub:

```bash
apptainer pull mrsa_proj.sif docker://nbisweden/workshop-reproducible-research
```

This should result in a SIF file called `mrsa_proj.sif`.

In the Docker image we included the code needed for the workflow in the
`/course` directory of the image. These files are of course also available in
the Apptainer image. However, a Apptainer image is read-only. This will be a
problem if we try to run the workflow within the `/course` directory, since the
workflow will produce files and Snakemake will create a `.snakemake` directory.
Instead, we need to provide the files externally from our host system and simply
use the Apptainer image as the environment to execute the workflow in (_i.e._
all the software and dependencies).

In your current working directory
(`workshop-reproducible-research/tutorials/containers/`) the vital MRSA project
files are already available (`Snakefile`, `config.yml` and
`code/supplementary_material.qmd`). Since Apptainer bind mounts the current
working directory we can simply execute the workflow and generate the output
files using:

```bash
apptainer run mrsa_proj.sif
```

This executes the default run command, which is `snakemake -rp -c 1 --configfile
config.yml` (as defined in the original `Dockerfile`). Once completed you should
see a bunch of directories and files generated in your current working
directory, including the `results/` directory containing the final HTML report.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to build a Apptainer image using Apptainer inside Docker.
> - How to convert Docker images to Apptainer images.
> - How to run Apptainer images.
