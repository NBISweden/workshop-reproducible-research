In the previous section we downloaded a Docker image of Ubuntu and noticed that
it was based on layers, each with a unique hash as id. An image in Docker is
based on a number of read-only layers, where each layer contains the
differences to the previous layers. If you've done the [Git tutorial](git-1-introduction)
this might remind you of how a Git commit contains the difference to the
previous commit. The great thing about this is that we can start from one base
layer, say containing an operating system and some utility programs, and then
generate many new images based on this, say 10 different project-specific
images. This dramatically reduces the storage space requirements. For example,
Bioconda (see the [Conda tutorial](conda-1-introduction)) has one base image
and then one individual layer for each of the more than 3000 packages
available in Bioconda.

Docker provides a convenient way to describe how to go from a base image to the
image we want by using a "Dockerfile". This is a simple text file containing
the instructions for how to generate each layer. Docker images are typically
quite large, often several GBs, while Dockerfiles are small and serve as
blueprints for the images. It is therefore good practice to have your
Dockerfile in your project Git repository, since it allows other users to
exactly replicate your project environment.

We will be looking at a Dockerfile called `Dockerfile_slim` that is located in your
`containers` directory (where you should hopefully be standing already). We will now
go through that file and discuss the different steps and what they do. After
that we'll build the image and test it out. Lastly, we'll start from that image
and make a new one to reproduce the results from the [Conda tutorial](conda-3-projects).

## Understanding Dockerfiles

Here are the first few lines of `Dockerfile_slim`. Each line in the Dockerfile
will typically result in one layer in the resulting image. The format for
Dockerfiles is `INSTRUCTION arguments`. A full specification of the format,
together with best practices, can be found
[here](https://docs.docker.com/engine/reference/builder/).

```Dockerfile
FROM condaforge/mambaforge

LABEL description = "Minimal image for the NBIS reproducible research course."
MAINTAINER "John Sundh" john.sundh@scilifelab.se
```

Here we use the instructions `FROM`, `LABEL` and `MAINTAINER`. While `LABEL`
and `MAINTAINER` is just meta-data that can be used for organizing your
various Docker components the important one is `FROM`, which specifies the
base image we want to start from. Because we want to use `mamba` to install
packages we will start from an image from the conda-forge community that has
`mamba` pre-installed. This image was in turn built using a Dockerfile as a
blueprint and then uploaded to
[Dockerhub](https://hub.docker.com/r/condaforge/mambaforge). The conda-forge
community keeps the Dockerfile in a git repository and you can view the file
[here](https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile).
You will see that it starts from an official Ubuntu image (check the first
line with the `FROM` instruction), followed by code to install various
packages including mamba.

> *Many roads to Rome* <br>
> When it comes to choosing the best image to start from there are multiple
> routes you could take. Say you want to run RStudio in a Conda environment
> through a Jupyter notebook. You could then start from one of the
> [rocker images](https://github.com/rocker-org/rocker) for R, a
> [Mambaforge image](https://hub.docker.com/r/condaforge/mambaforge), or
> a [Jupyter image](https://hub.docker.com/r/jupyter/). Or you just start
> from one of the low-level official images and set up everything from scratch.

Let's take a look at the next section of `Dockerfile_slim`.

```Dockerfile
# Use bash as shell
SHELL ["/bin/bash", "--login", "-c"]

# Set workdir
WORKDIR /course

# Set time zone
ENV TZ="Europe/Stockholm"
ENV DEBIAN_FRONTEND=noninteractive
```

`SHELL` simply sets which shell to use and `WORKDIR` determines the
directory the container should start in. The `ENV` instruction is used to
set environmental variables and here we use it to set the time zone by declaring
a `TZ` variable. The `DEBIAN_FRONTEND=noninteractive` line means that we
force the subsequent installation to not prompt us to set the time zone manually.

The next few lines introduce the important `RUN` instruction, which is used
for executing shell commands:

```Dockerfile
# Install package for setting time zone
RUN apt-get update && apt-get install -y tzdata && apt-get clean

# Configure Conda/Mamba
RUN mamba init bash && conda config --set channel_priority strict && \
    conda config --append channels bioconda && \
    conda config --append channels r && \
    conda config --set subdir linux-64
```

The first RUN command installs the `tzdata` package for managing local time
settings in the container. This may not always be required for your
Dockerfile but it's added here because some R packages used in the course
require it.

> **Note** <br>
> While installing things with `apt-get` inside Dockerfiles is relatively common
> practice, it's important to note that this _may_ affect reproducibility, since
> it's not common to specify an exact version. The packages installed in this
> manner are, however, usually not important for the actual analyses performed,
> but rather help in the building of the container image itself. While not
> critical, it's important to note this from a reproducibility perspective.

Next, we run `mamba init bash` to initialize the bash shell inside the
image, meaning we can use `mamba activate` in containers that run from the
image. In the same `RUN` statement we also configure the strict channel priority
and add appropriate channels with `conda config`. You'll probably recognize
this from the [pre-course-setup](../course-information/pre-course-setup).
The last part sets the somewhat obscure `subdir` config parameter pointing to
the `linux-64` architecture of conda channels.

As a general rule, you want each layer in an image to be a "logical unit". For
example, if you want to install a program the `RUN` command should both
retrieve the program, install it and perform any necessary clean up. This is
due to how layers work and how Docker decides what needs to be rerun between
builds. More on this later.

Next up is:

```Dockerfile
# Open port for running Jupyter Notebook
EXPOSE 8888

# Start Bash shell by default
CMD /bin/bash
```

`EXPOSE` opens up the port 8888, so that we can later run a Jupyter Notebook
server on that port. `CMD` is an interesting instruction. It sets what a
container should run when nothing else is specified, _i.e._ if you run
`docker run [OPTIONS] [IMAGE]` without the additional `[COMMAND] [ARG]`. It
can be used for example for printing some information on how to use the
image or, as here, start a Bash shell for the user. If the purpose of your
image is to accompany a publication then `CMD` could be to run the workflow that
generates the paper figures from raw data, _e.g._ `CMD snakemake -s
Snakefile -c 1 generate_figures`.

## Building from Dockerfiles

Now we understand how a Dockerfile works. Constructing the image itself from the
Dockerfile can be done as follows - try it out:

> **Important** <br>
> If your computer is a MAC with the M1 chip, you may have to add
> `--platform linux/x86_64` to the `docker build` command.

```bash
docker build -f Dockerfile_slim -t my_docker_image .
```

This should result in something similar to this:

```
 [+] Building 2.2s (7/7) FINISHED
 => [internal] load build definition from Dockerfile_slim                                                                                                                                             0.0s
 => => transferring dockerfile: 667B                                                                                                                                                                  0.0s
 => [internal] load .dockerignore                                                                                                                                                                     0.0s
 => => transferring context: 2B                                                                                                                                                                       0.0s
 => [internal] load metadata for docker.io/condaforge/mambaforge:latest                                                                                                                               0.0s
 => [1/3] FROM docker.io/condaforge/mambaforge                                                                                                                                                        0.0s
 => CACHED [2/3] WORKDIR /course                                                                                                                                                                      0.0s
 => [3/3] RUN mamba init bash && conda config --set channel_priority strict &&     conda config --append channels bioconda &&     conda config --append channels r &&     conda config --set subdir   2.1s
 => exporting to image                                                                                                                                                                                0.0s
 => => exporting layers                                                                                                                                                                               0.0s
 => => writing image sha256:53e6efeaa063eadf44c509c770d887af5e222151f08312e741aecc687e6e8981                                                                                                          0.0s
 => => naming to docker.io/library/my_docker_image
```

Exactly how the output looks depends on which version of Docker you are using.
The `-f` flag sets which Dockerfile to use and `-t` tags the image with a name.
This name is how you will refer to the image later. Lastly, the `.` is the path
to where the image should be build (`.` means the current directory). This had
no real impact in this case, but matters if you want to import files. Validate
with `docker image ls` that you can see your new image.

## Creating your own Dockerfile

Now it's time to make your own Dockerfile to reproduce the results from the
[Conda tutorial](conda-3-projects). If you haven't done the tutorial,
it boils down to creating a Conda environment file, setting up that
environment, downloading three RNA-seq data files, and running FastQC on those
files. We will later package and run the whole RNA-seq workflow in a Docker
container, but for now we keep it simple to reduce the size and time required.

The Conda tutorial uses a shell script, `run_qc.sh`, for downloading and
running the analysis. A copy of this file should also be available in your
current directory. If we want to use the same script we need to include it in
the image. A basic outline of what we need to do is:

1. Create a file called `Dockerfile_conda`
2. Start the image from the `my_docker_image` we just built
3. Install the package `fastqc` which is required for the analysis.
4. Add the `run_qc.sh` script to the image
5. Set the default command of the image to run the `run_qc.sh` script.

We'll now go through these steps in more detail. Try to add the
corresponding code to `Dockerfile_conda` on your own, and if you get stuck
you can click to reveal the solution below under "Click to show solution".

**Set image starting point**

To set the starting point of the new image, use the `FROM` instruction and
point to `my_docker_image` that we built in the previous _Building from
Dockerfiles_ step.

**Install packages**

Use the `RUN` instruction to install the package `fastqc=0.11.9` with Mamba.
Here there are several options available. For instance we could add an
environment file _e.g._ `environment.yml` from the Conda tutorial and use
`mamba env create` to create an environment from that file. Or we could
create an environment directly with `mamba create`. We'll try this later
option here, so add a line that will create an environment named
`project_mrsa` containing the two packages, and also clean up packages and
cache after installation. Use the `-y` flag to `mamba create` to avoid the
prompt that expects an interaction from the user.

In order to have the `project_mrsa` environment activated upon start-up we
need to add two more lines to the Dockerfile. First we need to use a `RUN`
instruction to run `echo "source activate project_mrsa" >> ~/.bashrc`, and
then we need to use the `ENV` instruction to set the `$PATH` variable
inside the image to `/opt/conda/envs/project_mrsa/bin:$PATH`.

**Add the analysis script**

Use the `COPY` instruction to Add `run_qc.sh` to the image. The syntax is
`COPY SOURCE TARGET`. In this case `SOURCE` is the `run_qc.sh`
script and `TARGET` is a path inside the image, for simplicity it can be
specified with `./`.

**Set default command**

Use the `CMD` instruction to set the default command for the image to `bash
run_qc.sh`.

<details>
<summary> Click to show </summary>

```Dockerfile
FROM my_docker_image

RUN mamba create -y -n project_mrsa -c bioconda fastqc=0.11.9 && mamba clean -a

RUN echo "source activate project_mrsa" >> ~/.bashrc

ENV PATH=/opt/conda/envs/project_mrsa/bin:$PATH

COPY run_qc.sh .

CMD bash run_qc.sh
```

</details>

Build the image and tag it `my_docker_conda` (remember to add
`--platform linux/x86_64` to the build command if you are using a Mac with the
Apple chip).

```bash
docker build -t my_docker_conda -f Dockerfile_conda .
```

Verify that the image was built using `docker image ls`.

> **Quick recap** <br>
> In this section we've learned:
>
> - How the keywords `FROM`, `LABEL`, `MAINTAINER`, `RUN`, `ENV`, `SHELL`,
>  `WORKDIR`, and `CMD` can be used when writing a Dockerfile.
> - How to use `docker build` to construct and tag an image from a Dockerfile.
> - How to create your own Dockerfile.
