# Building images

Singularity has the ability to convert Docker images to the Singularity Image
Format (SIF). We can try this out by running:

```bash
singularity pull docker://godlovedc/lolcow
```

This command generates a .sif file where the individual layers of the specified
Docker image have been combined and converted to Singularity's native format.
We can now use `run`, `exec`, and `shell` commands on this image file. Try it:

```bash
singularity run lolcow_latest.sif
singularity exec lolcow_latest.sif fortune
singularity shell lolcow_latest.sif
```

> **Quick recap** <br>
> In this section we covered how to use `singularity pull` to download and
> run Docker images as Singularity containers.

## Building from scratch

As we have seen, it is possible to convert Docker images to the Singularity
format when needed and run them using Singularity. In terms of making
a research project reproducible using containers, it may be enough to *e.g.*
define a Dockerfile (recipe for a Docker image) as well as supply a Docker
image for others to download and use, either directly through Docker, or by
Singularity. Even better, from a reproducibility aspect, would be to also
generate the Singularity image from the Docker image and provide that for
potential future users (since the image is a static file, whereas running
`singularity pull` or `singularity build` would rebuild the image at the time
of issuing the command).

A third option, would be to define a Singularity recipe, either on its own or
in addition to the Dockerfile. The equivalent of a Dockerfile for Singularity
is called a Singularity Definition file ("def file").

The def file consists of two parts:

* A header that defines the core OS and related features
* Optional sections, each starting with a `%`-sign, that add content or execute
  commands during the build process

As an example, we can look at the def file used for the image we played with
above in the `singularity/` directory (we previously pulled lolcow from
Dockerhub, but it also exists in the Singularity library and can be pulled by
`singularity pull library://godlovedc/funny/lolcow`). The lolcow def file looks
like this:

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

The first part of the header sets the bootstrap agent. In the lolcow example
DockerHub is used. Alternatively one could set it to *library* to use the
Singularity Library. There are also other bootstrap agents available (see [this
link](https://sylabs.io/guides/3.3/user-guide/definition_files.html#preferred-bootstrap-agents)
for details). Next, the base image that the new image starts from is defined,
in this case the Docker image `ubuntu:16.04`.

In the lolcow def file three sections are used (`%post`, `%environment`, and
`runscript`).

* `%post` is similar to the `RUN` instruction in Dockerfiles. Here is where you
  include code to download files from the internet, install software, create
  directories etc.
* `%environment` is similar to the `ENV` instruction in Dockerfiles. It is used
  to set environmental variables that will be available when running the
  container. The variables set in this section will not however be available
  during the build and should in the cases they are needed then also be set in
  the `%post` section.
* `%runscript` is similar to the `CMD` instruction in Dockerfiles and contains
  the default command to be executed when running the container.

> **Quick recap** <br>
> In this section we covered the basic parts of a Singularity definition
> file (a def file), including  `BootStrap`, `From` `%post`, `%environment`
> and `%runscript`.

## Singularity def file for the MRSA project

Let's use the MRSA case study project to define our own Singularity def file!
We will not make an image for the whole workflow but rather focus on the
`run_qc.sh` script that we used in the end of the [conda tutorial](conda.md).
This script is included in the `code` directory of your current working
directory (`singularity`) and, when executed, downloads a few fastq-files and
runs FastQC on them. To run the script we need the software SRA-Tools and
FastQC.

* Make a new file called `run_qc.def` and add the following lines:

```
Bootstrap: library
From: ubuntu:16.04

%labels
    DESCRIPTION Image for running the run_qc.sh script
    AUTHOR <Your Name>
```

Here we'll use the Singularity Library as bootstrap agent, instead of DockerHub 
as in the lol_cow example above. The base Singularity image will be 
`ubuntu:16.04`. We can also add metadata to our image using any name-value pair.

* Next, add the `%environment` section:

```
%environment
    export LC_ALL=C
    export PATH=/usr/miniconda3/bin:$PATH
```

This sets the default locale as well as includes the PATH to conda (which we
will soon install).

* Now add the `%post` section:

```
%post
    apt-get update
    apt-get install -y --no-install-recommends bzip2 ca-certificates curl
    apt-get clean

    # Install conda:
    curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O
    bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/
    rm Miniconda3-4.7.12.1-Linux-x86_64.sh

    # Configure conda:
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # Install requirements:
    conda install -c bioconda fastqc=0.11.9 sra-tools=2.10.0
    conda clean --all
```

You should recognize parts of this from the Docker tutorial. Basically, we 
install some required basic tools like bzip2, ca-certificates and curl, then 
install and configure conda and finally install the required tools for the 
`run_qc.sh` script.

* Next, add a `%test` section:

```
%test
    fastqc --version
    fastq-dump --version
```

The test section runs at the end of the build process and you can include any
code here to verify that your image works as intended. Here we just make sure
that the `fastqc` and `fastq-dump` are available.

* Finally, add the `%runscript`:

```
%runscript
    bash code/run_qc.sh
```

We should now be ready to build our image from this def file using 
`singularity build`. Now, depending on the system you are running on and the 
version of Singularity, you may not have the option to build locally. However, 
Singularity has the option to build images remotely. To do this, you need to:

* Go to [https://cloud.sylabs.io/library](https://cloud.sylabs.io/library) and 
  create an account
* Log in and find "Access Tokens" in the menu and create a new token
* Copy the token
* In your terminal, run `singularity remote login` and hit ENTER. You should be
  asked to enter the token (API Key). Paste the copied token and hit ENTER. 
  You should get a **API Key Verified!** message.

> **Attention** <br>
> In case you are not asked to enter the API Key, you can try to run 
> `singularity remote login SylabsCloud` instead.

We can now try to build the MRSA Singularity image using the `--remote` flag:

```bash
singularity build --remote run_qc.sif run_qc.def
```

Did it work? Can you figure out why it failed? Tip: it has to do with the PATH
(well, to be fair, it could fail for several reasons, but if you did everything
else correct it should be related to the PATH).

??? note "Click to see the solution"
    You need to add conda to the PATH. `%environment` makes it available at
    runtime but not during build.

    Update the `%post` section as follows:

    ```bash
    # Install conda:
    curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O
    bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/
    rm Miniconda3-4.7.12.1-Linux-x86_64.sh
    export PATH=/usr/miniconda3/bin:$PATH ## <- add this line
    ```

    You also need to update the `%test` section:

    ```bash
    export PATH=/usr/miniconda3/bin:$PATH ## <- add this line
    fastqc --version
    fastq-dump --version
    ```

The build should now hopefully work and produce a Singularity image called
`run_qc.sif`. To run the image, *i.e.* executing `code/run_qc.sh` using the
tools in the container, do the following:

```bash
singularity run run_qc.sif
```

The fastq-files should now be downloaded and FastQC should be run on these
files producing output directories and files in your current working directory.

> **Tip** <br>
> We do not cover it here but it is possible to build Singularity images as
> writable sandbox images. This enables starting a shell in the container and
> *e.g.* installing software. This may be convenient during the design of the
> definition file to test what commands to include. When everything is
> working as expected one should rebuild directly from the definition file to
> a final SIF file.

> **Note** <br>
> A somewhat important section that we have not covered here is the `%files`
> section. This is similar to the `ADD` or `COPY` instructions in
> a Dockerfile. One simply defines, on each line, a file to be copied from
> host to the container image using the format `<source> <destination>`. This
> does not currently work with `--remote` building though.
