We're almost ready to start, just one last note on nomenclature. You might have
noticed that we sometimes refer to "Docker images" and sometimes to "Docker
containers". We use images to start containers, so containers are simply an
instances of an image. You can have an image containing, say, a certain
Linux distribution, and then start multiple containers running that
same OS.

> **Attention!** <br>
> If you don't have root privileges you have to prepend all Docker commands
> with `sudo`.

## Downloading images

Docker containers typically run Linux, so let's start by downloading an image
containing Ubuntu (a popular Linux distribution that is based on only
open-source tools) through the command line.

```bash
docker pull ubuntu:latest
```

You will notice that it downloads different layers with weird hashes as names.
This represents a very fundamental property of Docker images that we'll get
back to in just a little while. The process should end with something along the
lines of:

```no-highlight
Status: Downloaded newer image for ubuntu:latest
docker.io/library/ubuntu:latest
```

Let's take a look at our new and growing collection of Docker images:

```bash
docker image ls
```

The Ubuntu image should show up in this list, with something looking like this:

```
REPOSITORY       TAG              IMAGE ID            CREATED             SIZE
ubuntu           latest           d70eaf7277ea        3 weeks ago         72.9MB
```

## Running containers

We can now start a container from the image we just downloaded. We can refer to
the image either by "REPOSITORY:TAG" ("latest" is the default so we can omit
it) or "IMAGE ID". The syntax for `docker run` is `docker run [OPTIONS] IMAGE
[COMMAND] [ARG...]`. To see the available options run `docker run --help`.
The `COMMAND` part is any command that you want to run inside the container,
it can be a script that you have written yourself, a command line tool or a
complete workflow. The `ARG` part is where you put optional arguments that
the command will use.

Let's run `uname -a` to get some info about the operating system. In this
case, `uname` is the `COMMAND` and `-a` the `ARG`. This command will display
some general info about your system, and the `-a` argument tells `uname` to
display all possible information.

First run it on your own system (use `systeminfo` if you are on Windows):

```bash
uname -a
```

This should print something like this to your command line:

```no-highlight
Darwin liv433l.lan 15.6.0 Darwin Kernel Version 15.6.0: Mon Oct  2 22:20:08 PDT 2017; root:xnu-3248.71.4~1/RELEASE_X86_64 x86_64
```

Seems like I'm running the Darwin version of MacOS. Then run it in the Ubuntu
Docker container:

```bash
docker run ubuntu uname -a
```

Here I get the following result:

```no-highlight
Linux 24d063b5d877 5.4.39-linuxkit #1 SMP Fri May 8 23:03:06 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
```

And now I'm running on Linux! What happens is that we use the downloaded
ubuntu image to run a container that has `Ubuntu` as the operating system,
and we instruct Docker to execute `uname -a` to print the system info within
that container. The output from the command is printed to the terminal.

Try the same thing with `whoami`.

## Running interactively

So, seems we can execute arbitrary commands on Linux. This looks useful, but
maybe a bit limited. We can also get an interactive terminal with the flags
`-it`.

```bash
docker run -it ubuntu
```

Your prompt should now look similar to:

```no-highlight
root@1f339e929fa9:/#
```

You are now using a terminal inside a container running Ubuntu. Here you can
do whatever; install, run, remove stuff. Anything you do will be isolated
within the container and never affect your host system.

Now exit the container with `exit`.

## Containers inside scripts

Okay, so Docker lets us work in any OS in a quite convenient way. That would
probably be useful on its own, but Docker is much more powerful than that. For
example, let's look at the `shell` part of the `index_genome` rule in
the Snakemake workflow for the MRSA case study:

```python
shell:
    """
    bowtie2-build tempfile intermediate/{wildcards.genome_id} > {log}
    """
```

You may have seen that one can use containers through both Snakemake and
Nextflow if you've gone through their tutorial's extra material, but we can
also use containers directly inside scripts in a very simple way. Let's imagine
we want to run the above command using containers instead. How would that look?
It's quite simple, really: first we find a container image that has `bowtie2`
installed, and then prepend the command with `docker run <image>`.

First of all we need to download the genome to index though, so run:
```bash
curl -o NCTC8325.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
gunzip -c NCTC8325.fa.gz > tempfile
```

To download and prepare the input for Bowtie2.

Now try running the following Bash code:

```bash
docker run -v $(pwd):/analysis quay.io/biocontainers/bowtie2:2.5.1--py39h3321a2d_0 bowtie2-build /analysis/tempfile /analysis/NCTC8325
```

Docker will automatically download the container image for Bowtie2 version
2.5.1 from the remote repository `https://quay.io/repository/biocontainers/bowtie2`
and subsequently run the command! This is the `docker run [OPTIONS] IMAGE
[COMMAND] [ARG...]` syntax just like before. In this case
`quay.io/biocontainers/bowtie2:2.5.1--py39h3321a2d_0` is the IMAGE but
instead of first downloading and then running it we point to its remote
location directly, which will cause Docker to download it on
the fly. The `bowtie2-build` part is the COMMAND followed by the ARG (the
input tempfile and the output index)

The `-v $(pwd):/analysis` part is the OPTIONS which we use to mount the current
directory inside the container in order to make the `tempfile` input available
to Bowtie2. More on these so-called "Bind mounts" in Section 4 of this tutorial.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to use `docker pull` for downloading remotely stored images
> - How to use `docker image ls` for getting information about the images we
>   have on our system.
> - How to use `docker run` for starting a container from an image.
> - How to use the `-it` flag for running in interactive mode.
> - How to use Docker inside scripts.
