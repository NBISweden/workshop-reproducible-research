<script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML'></script>

In the previous section we downloaded a Docker image of Ubuntu and noticed that
it was based on layers, each with a unique hash as id. An image in Docker is
based on a number of read-only layers, where each layer contains the
differences to the previous layers. If you've done the [Git tutorial](git-1-introduction)
this might remind you of how a Git commit contains the difference to the
previous commit. The great thing about this is that we can start from one base
layer, say containing an operating system and some utility programs, and then
generate many new images based on this, say 10 different project-specific
images. The total space requirements would then only be $base+\sum_{i=1}^{10}(specific_{i})$
rather than $\sum_{i=1}^{10}(base+specific_{i})$. For example, Bioconda (see the
[Conda tutorial](conda-1-introduction)) has one base image and then
one individual layer for each of the more than 3000 packages available in
Bioconda.

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

```no-highlight
FROM ubuntu:16.04

LABEL description = "Minimal image for the NBIS reproducible research course."
MAINTAINER "John Sundh" john.sundh@scilifelab.se
```

Here we use the instructions `FROM`, `LABEL` and `MAINTAINER`. The important
one is `FROM`, which specifies the base image our image should start from. In
this case we want it to be Ubuntu 16.04, which is one of the [official
repositories](https://hub.docker.com/explore/). There are many roads to Rome
when it comes to choosing the best image to start from. Say you want to run
RStudio in a Conda environment through a Jupyter notebook. You could then start
from one of the [rocker images](https://github.com/rocker-org/rocker) for R,
a [Miniconda image](https://hub.docker.com/r/continuumio/miniconda/), or
a [Jupyter image](https://hub.docker.com/r/jupyter/). Or you just start from
one of the low-level official images and set up everything from scratch.
`LABEL` and `MAINTAINER` is just meta-data that can be used for organizing your
various Docker components.

Let's take a look at the next section of the Dockerfile.

```no-highlight
# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /course
```

`SHELL` simply sets which shell to use. `WORKDIR` determines the directory the
container should start in. The next few lines introduce the important `RUN`
instruction, which is used for executing shell commands:

```no-highlight
# Install necessary tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends bzip2 \
                                               ca-certificates \
                                               curl \
                                               fontconfig \
                                               git \
                                               language-pack-en \
                                               tzdata \
                                               vim \
                                               unzip \
                                               wget \
    && apt-get clean

# Install Miniconda and add to PATH
RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O && \
    bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/ && \
    rm Miniconda3-4.7.12.1-Linux-x86_64.sh && \
    /usr/miniconda3/bin/conda clean -tipsy && \
    ln -s /usr/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /usr/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
```

As a general rule, you want each layer in an image to be a "logical unit". For
example, if you want to install a program the `RUN` command should both
retrieve the program, install it and perform any necessary clean up. This is
due to how layers work and how Docker decides what needs to be rerun between
builds. The first command uses Ubuntu's package manager APT to install some
packages (similar to how we've previously used Conda). Say that the first
command was split into two instead:

```no-highlight
# Update apt-get
RUN apt-get update

# Install packages
RUN apt-get install -y --no-install-recommends bzip2 \
                                               ca-certificates \
                                               curl \
                                               fontconfig \
                                               git \
                                               language-pack-en \
                                               tzdata \
                                               vim \
                                               unzip \
                                               wget

# Clear the local repository of retrieved package files
RUN apt-get clean
```

The first command will update the apt-get package lists and the second will
install the packages `bzip2`, `ca-certificates`, `curl`, `fontconfig`, `git`,
`language-pack-en`, `tzdata`, `vim`, `unzip` and `wget`. Say that you build this
image now, and then in a month's time you realize that you would have liked a
Swedish language pack instead of an English. You change to `language-pack-sv`
and rebuild the image. Docker detects that there is no layer with the new
list of packages and reruns the second `RUN` command. *However, there is no
way for Docker to know that it should also update the apt-get package lists*.
You therefore risk to end up with old versions of packages and, even worse,
the versions would depend on when the previous version of the image was first
built.

The next `RUN` command retrieves and installs Miniconda3. Let's see what would
happen if we had that as separate commands instead.

```no-highlight
# Download Miniconda3
RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O

# Install it
RUN bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/

# Remove the downloaded installation file
RUN rm Miniconda3-4.7.12.1-Linux-x86_64.sh

# Remove unused packages and caches
RUN /usr/miniconda3/bin/conda clean -tipsy

# Permanently enable the Conda command
RUN ln -s /usr/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh
RUN echo ". /usr/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc

# Add the base environment permanently to PATH
RUN echo "conda activate base" >> ~/.bashrc
```

Remember that each layer contains the difference compared to the previous
layer? What will happen here is that the first command adds the installation
file and the second will unpack the file and install the software. The third
layer will say "the installation file should no longer exist on the file
system". However, the file will still remain in the image since the image is
constructed layer-by-layer bottom-up. This results in unnecessarily many layers
and bloated images. Line four is cleaning up conda to free up space, and the
next two lines are there to make the Conda command available in the shell.
The last command adds a code snippet to the bash startup file which
automatically activates the Conda base environment in the container.

```no-highlight
# Add conda to PATH and set locale
ENV PATH="/usr/miniconda3/bin:${PATH}"
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
```

Here we use the new instruction `ENV`. The first command adds `conda` to the
path, so we can write `conda install` instead of `/usr/miniconda3/bin/conda install`.
The next two commands set an UTF-8 character encoding so that we can use
weird characters (and a bunch of other things).

```no-highlight
# Configure Conda channels and install Mamba
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda install mamba \
    && mamba clean --all
```

Here we just configure Conda and install Mamba, for quicker installations of
any subsequent Conda packages we might want to do.

```no-highlight
# Open port for running Jupyter Notebook
EXPOSE 8888

# Start Bash shell by default
CMD /bin/bash
```

`EXPOSE` opens up the port 8888, so that we can later run a Jupyter Notebook
server on that port. `CMD` is an interesting instruction. It sets what a
container should run when nothing else is specified. It can be used for example
for printing some information on how to use the image or, as here, start a shell
for the user. If the purpose of your image is to accompany a publication then
`CMD` could be to run the workflow that generates the paper figures from raw
data.

## Building from Dockerfiles

Ok, so now we understand how a Dockerfile works. Constructing the image from
the Dockerfile is really simple. Try it out now:

```bash
docker build -f Dockerfile_slim -t my_docker_image .
```

This should result in something similar to this:

```
=> [internal] load build definition from Dockerfile_slim
=> => transferring dockerfile: 1.88kB
=> [internal] load .dockerignore
=> => transferring context: 2B
=> [internal] load metadata for docker.io/library/ubuntu:16.04
=> [auth] library/ubuntu:pull token for registry-1.docker.io
=> [1/5] FROM docker.io/library/ubuntu:16.04@sha256:bb84bbf2ff36d46acaf0bb0c6bcb33dae64cd93cba8652d74c9aaf438fada438
=> CACHED [2/5] WORKDIR /course
=> CACHED [3/5] RUN apt-get update &&     apt-get install -y --no-install-recommends bzip2                                                ca-certificates
=> CACHED [4/5] RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O &&   h Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3
=> CACHED [5/5] RUN conda config --add channels bioconda     && conda config --add channels conda-forge     && conda config --set channel_priority strict     && conda instal
=> exporting to image
=> => exporting layers
=> => writing image sha256:d14301f829d4554816df54ace927ec0aaad4a994e028371455f7a18a370f6af9
=> => naming to docker.io/library/my_docker_image
```

Exactly how the output looks depends on which version of Docker you are using.
The `-f` flag sets which Dockerfile to use and `-t` tags the image with a name.
This name is how you will refer to the image later. Lastly, the `.` is the path
to where the image should be build (`.` means the current directory). This had
no real impact in this case, but matters if you want to import files. Validate
with `docker image ls` that you can see your new image.

## Creating your own Dockerfile

Now it's time to make our own Dockerfile to reproduce the results from the
[Conda tutorial](conda-3-projects). If you haven't done the tutorial,
it boils down to creating a Conda environment file, setting up that
environment, downloading three RNA-seq data files, and running FastQC on those
files. We will later package and run the whole RNA-seq workflow in a Docker
container, but for now we keep it simple to reduce the size and time required.

The Conda tutorial uses a shell script, `run_qc.sh`, for downloading and
running the analysis. A copy of this file should also be available in your
current directory. If we want to use the same script we need to include it in
the image. So, this is what we need to do:

1. Create the file `Dockerfile_conda`.

2. Set `FROM` to the image we just built.

3. Install the required packages with Conda. We could do this by adding
   `environment.yml` from the Conda tutorial, but here we do it directly as
   `RUN` commands. We need to add the conda-forge and bioconda channels with
   `conda config --add channels <channel_name>` and install `fastqc=0.11.9` and
   `sra-tools=2.10.1` with `conda install`. The packages will be installed to
   the default environment named `base` inside the container.

4. Add `run_qc.sh` to the image by using the `COPY` instruction. The syntax is
   `COPY source target`, so in our case simply `COPY run_qc.sh .` to copy to
   the work directory in the image.

5. Set the default command for the image to `bash run_qc.sh`, which will
   execute the shell script.

Try to add required lines to `Dockerfile_conda`. If it seems overwhelming you
can take a look at an example below:

??? example "Click to show the solution"
    ```no-highlight
    FROM my_docker_image:latest
    RUN conda config --add channels bioconda && \
        conda config --add channels conda-forge && \
        mamba install -n base fastqc=0.11.9 sra-tools=2.10.1
    COPY run_qc.sh .
    CMD bash run_qc.sh
    ```

Build the image and tag it `my_docker_conda`:

```bash
docker build -t my_docker_conda -f Dockerfile_conda .
```

Verify that the image was built using `docker image ls`.

!!! Success "Quick recap"
    In this section we've learned:

    - How the keywords `FROM`, `LABEL`, `MAINTAINER`, `RUN`, `ENV`, `SHELL`,
    `WORKDIR`, and `CMD` can be used when writing a Dockerfile.
    - The importance of letting each layer in the Dockerfile be a "logical unit".
    - How to use `docker build` to construct and tag an image from a Dockerfile.
    - How to create your own Dockerfile.
