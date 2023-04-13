## Singularity as an alternative container tool

Singularity is a container software alternative to Docker. It was originally
developed by researchers at Lawrence Berkeley National Laboratory with focus on
security, scientific software, and HPC clusters. One of the ways in which
Singularity is more suitable for HPC is that it very actively restricts
permissions so that you do not gain access to additional resources while inside
the container. Singularity also, unlike Docker, stores images as single files
using the *Singularity Image Format* (SIF). A SIF file is self-contained and can
be moved around and shared like any other file, which also makes it easy to work
with on an HPC cluster.

> **Singularity and Apptainer** <br>
> The open source Singularity project was recently renamed to *Apptainer*. 
> Confusingly, the company *Sylabs* still keeps their commercial branch of 
> the project under the Singularity name, and offer a free 'Community 
> Edition' version. The name change was done in order to clarify the 
> distinction between the open source project and the various commercial 
> versions.
> At the moment there is virtually no difference to you as a user whether you 
> use Singularity or Apptainer, but eventually it's very likely that the two 
> will diverge. 
> We have opted to stick with the original name in the material for now, 
> while the change is still being adopted by the community and various 
> documentation online. In the future we will however move to using only 
> Apptainer to follow the open source route of the project.

While it is possible to define and build Singularity images from scratch, in a
manner similar to what you've already learned for Docker, this is not something
we will cover here (but feel free to read more about this in *e.g.* the
[Singularity docs](https://sylabs.io/guides/master/user-guide/) or the [Uppmax
Singularity user guide](https://www.uppmax.uu.se/support/user-guides/singularity-user-guide/)).

The reasons for not covering Singularity more in-depth are varied, but it
basically boils down to it being more or less Linux-only, unless you use Virtual
Machines (VMs). Even with this you'll run into issues of incompatibility of
various kinds, and these issues are further compounded if you're on one of the 
new ARM64-Macs. You also need `root` (admin) access in order to actually *build*
Singularity images regardless of platform, meaning that you can't build them on
*e.g.* Uppmax, even though Singularity is already installed there. You can,
however, use the `--remote` flag, which runs the build on Singularity's own
servers. This doesn't work in practice a lot of the time, though, since most
scientist will work in private Git repositories so that their research and code
is not available to anybody, and the `--remote` flag requires that *e.g.* the
`environment.yml` file is publicly available.

There are very good reasons to use Singularity, however, the major one being
that you aren't allowed to use Docker on most HPC systems! One of the nicer
features of Singularity is that it can convert Docker images directly for use
within Singularity, which is highly useful for the cases when you already built
your Docker image or if you're using a remotely available image stored on *e.g.*
DockerHub. For a lot of scientific work based in R and/or Python, however, it is
most often the case that you build your own images, since you have a complex
dependency tree of software packages not readily available in existing images.
So, we now have another problem for building our own images:

1. Only Singularity is allowed on HPC systems, but you can't build images there
   due to not having `root` access.
2. You can build Singularity images locally and transfer them to HPCs, but this
   is problematic unless you're running Linux natively.

Seems like a "catch 22"-problem, right? There are certainly workarounds (some of
which we have already mentioned) but most are roundabout or difficult to get
working for all use-cases. Funnily enough, there's a simple solution: run
Singularity locally from inside a Docker container! Conceptually very meta, yes,
but works very well in practice. What we are basically advocating for is that
you stick with Docker for most of your container-based work, but convert your
Docker images using Singularity-in-Docker whenever you need to work on an HPC.
This is of course not applicable to Linux users or those of you who are fine
with working through using VMs and managing any issues that arise from doing
that.

> **Summary** <br>
> Singularity/Apptainer is a great piece of software that is easiest to use if
> you're working on a Linux environment. Docker is, however, easier to use from
> a cross-platform standpoint and covers all use-cases except running on HPCs.
> Running on HPCs can be done by converting existing Docker images at runtime,
> while building images for use on HPCs can be done using local Docker images
> and Singularity-in-Docker.

## Singularity-in-Docker

By creating a bare-bones, Linux-based Docker image with Singularity you can
build Singularity images locally on non-Linux operating systems. 

> **Mac M-chip users** <br>
> If you're using a Mac laptop with the new M-chips (M1 or M2) then scroll 
> down to the **Singularity for ARM64 architecture** section.  

There is already a good image setup for just this, and it is defined in this 
[GitHub repository](https://github.com/kaczmarj/apptainer-in-docker). 
Looking at the instructions there we can see that we need to do the following:

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
while the last line is the Singularity/Apptainer command that actually does 
the conversion. All we need to do is to replace the `<IMAGE>` part with the 
Docker image we want to convert, *e.g.* `my_docker_image`.

* Replace `<IMAGE>` and `<TAG>` with one of your locally available Docker images
  and one of its tags and run the command - remember that you can use `docker
  image ls` to check what images you have available.

In the end you'll have a SIF file (*e.g.* `my_docker_image.sif`) that you can
transfer to an HPC such as Uppmax and run whatever analyses you need. If you
want to be able to do this without having to remember all the code you can check
out the [this script](https://github.com/fasterius/dotfiles/blob/main/scripts/singularity-in-docker.sh).

### Singularity for ARM64 architecture

The new M-chip series of Mac laptops have a CPU based on the 
[ARM](https://en.wikipedia.org/wiki/ARM_architecture_family) architecture 
which differs from the 
[AMD](https://en.wikipedia.org/wiki/List_of_AMD_CPU_microarchitectures) 
architecture. There's enough difference between these two systems that most
programs built for one will not work on the other. This causes issues when 
using images on one system that is built for the other system, and why you 
may have had to use the `--platform linux/x86_64` flag in the 
[Building images](containers-3-building-images.md) section of this tutorial. 
Docker is pretty good at handling this so that you can build an image on one 
system (_e.g._ ARM) for use on another system (_e.g._ AMD). However, at the 
time of writing this tutorial it is not possible to both build **and run** a 
Singularity/Apptainer image built for the AMD architecture on a system with the 
ARM architecture. To get around this we have designed a Docker image with 
Singularity built specifically for ARM (M-chip macs) and also created a 
workflow that uses packages that we know work on the ARM architecture.

So if you're on one of the new M-chip Macs, 


## Running Singularity

The following exercises assume that you have a login to the Uppmax HPC cluster
in Uppsala, but will also work for any other system that has Singularity
installed - like if you managed to install Singularity on your local system or
have access to some other HPC cluster. Let's try to convert the Docker image for
this course directly from DockerHub:

```bash
singularity pull mrsa_proj.sif docker://nbisweden/workshop-reproducible-research
```

This should result in a SIF file called `mrsa_proj.sif`.

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
`code/header.tex` and `code/supplementary_material.Rmd`). Since Singularity bind
mounts the current working directory we can simply execute the workflow and
generate the output files using:

```bash
singularity run mrsa_proj.sif
```

This executes the default run command, which is `snakemake -rp -c 1 --configfile
config.yml` (as defined in the original `Dockerfile`). Once completed you should
see a bunch of directories and files generated in your current working
directory, including the `results/` directory containing the final HTML report.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to build a Singularity image using Singularity inside Docker.
> - How to convert Docker images to Singularity images.
> - How to run Singularity images.
