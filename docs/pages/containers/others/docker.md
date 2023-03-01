# Introduction

Docker is a tool designed to make it easier to create, deploy, and run
applications by isolating them in "containers". The idea is to package your
code together with everything it needs (other packages it depends on, various
environment settings, data..) into one unit, *i.e.* a container. This way we
can ensure that the code results in exactly the same output regardless of where
it's executed. Containers are in many ways similar to virtual machines but more
lightweight. Rather than starting up a whole new OS, Docker containers can use
the same Linux kernel as the system that they're running on. This makes them
much faster and smaller compared to virtual machines. While this might sound
a bit technical, actually using Docker is quite easy, fun and very powerful.

Just as with Git, Docker was designed for software development but is rapidly
becoming widely used in scientific research. Say that you are building a web
application. You could then run the web server in one container and the
database in another, thereby reducing the risk of one system affecting the
other in unpredictable ways. Docker containers have also proven to be a very
good solution for packaging, running and distributing scientific data analyses.
Some applications relevant for reproducible research can be:

* When publishing, package your whole analysis pipeline together with your data
  in a Docker image and let it accompany the article. This way interested
  readers can reproduce your analysis at the push of a button.
* Packaging your analysis in a Docker container enables you to develop on
  *e.g.* your laptop and then seamlessly move to cluster or cloud to run the
  actual analysis.
* Say that you are collaborating on a project and you are using Mac while your
  collaborator is using Windows. You can then set up a Docker image specific
  for your project to ensure that you are working in an identical environment.

All of this might sound a bit abstract so far, but it'll become more clear
during the exercises. If you want to read more you can check out these
resources:

* A "Get started with Docker" at
  [docker.com](https://docs.docker.com/get-started/).
* An [early paper](https://arxiv.org/abs/1410.0846) on the subject of using
  Docker for reproducible research.

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](setup.md) for instructions on how to set it up if you haven't done so
already. Then open up a terminal and go to `workshop-reproducible-research/docker`.

!!! attention
    Docker images tend to take up quite a lot of space. In order to do all the
    exercises in this tutorial you need to have ~10 GB available.

## The basics

We're almost ready to start, just one last note on nomenclature. You might have
noticed that we sometimes refer to "Docker images" and sometimes to "Docker
containers". A container is simply an instance of an image. To use
a programming metaphor, if an image is a class, then a container is an instance
of that class — a runtime object. You can have an image containing, say,
a certain Linux distribution, and then start multiple containers running that
same OS.

!!! attention
    If you don't have root privileges you have to prepend all Docker commands
    with `sudo ` (see the tip above)

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
docker images
```

The Ubuntu image show show up in this list, with something looking like this:

```
REPOSITORY       TAG              IMAGE ID            CREATED             SIZE
ubuntu           latest           d70eaf7277ea        3 weeks ago         72.9MB
```

We can now start a container running our image. We can refer to the image
either by "REPOSITORY:TAG" ("latest" is the default so we can omit it) or
"IMAGE ID". The syntax for `docker run` is `docker run [OPTIONS] IMAGE
[COMMAND] [ARG...]`. Let's run the command `uname -a` to get some info about
the operating system. First run on your own system (skip this if you're using
Windows via the Windows 10 PowerShell, or use `systeminfo` which is the
Windows equivalent):

```bash
uname -a
```

This should print something like this to your command line:

```no-highlight
Darwin liv433l.lan 15.6.0 Darwin Kernel Version 15.6.0: Mon Oct  2 22:20:08 PDT 2017; root:xnu-3248.71.4~1/RELEASE_X86_64 x86_64
```

Seems like I'm running the Darwin version of macOS. Then run it in the Ubuntu
Docker container:

```bash
docker run ubuntu uname -a
```

Here I get the following result:

```no-highlight
Linux 24d063b5d877 5.4.39-linuxkit #1 SMP Fri May 8 23:03:06 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
```

And now I'm running on Linux! Try the same thing with `whoami`.

So, seems we can execute arbitrary commands on Linux. Seems useful, but maybe
a bit limited. We can also get an interactive terminal with the flags `-it`.

```bash
docker run -it ubuntu
```

This should put at a terminal prompt inside a container running Ubuntu. Your
prompt should now look similar to:

```no-highlight
root@1f339e929fa9:/#
```

Here you can do whatever; install, run, remove stuff. It will still be within
the container and never affect your host system. Now exit the container with
`exit`.

Ok, so Docker lets us work in any OS in a quite convenient way. That would
probably be useful on its own, but Docker is much more powerful than that.

!!! note "Quick recap"
    In this section we've learned:

    * How to use `docker pull` for downloading images from a central registry.
    * How to use `docker images` for getting information about the images we
      have on our system.
    * How to use `docker run` for starting a container from an image.
    * How to use the `-it` flag for running in interactive mode.

## Building a Docker image

<script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML'></script>
In the previous section we downloaded a Docker image of Ubuntu and noticed that
it was based on layers, each with a unique hash as id. An image in Docker is
based on a number of read-only layers, where each layer contains the
differences to the previous layers. If you've done the [Git tutorial](git.md)
this might remind you of how a Git commit contains the difference to the
previous commit. The great thing about this is that we can start from one base
layer, say containing an operating system and some utility programs, and then
generate many new images based on this, say 10 different project-specific
images. The total space requirements would then only be
$base+\sum_{i=1}^{10}(specific_{i})$ rather than
$\sum_{i=1}^{10}(base+specific_{i})$. For example, Bioconda (see the [Conda
tutorial](conda.md)) has one base image and then one individual layer for each
of the >3000 packages available in Bioconda.

Docker provides a convenient way to describe how to go from a base image to the
image we want by using a "Dockerfile". This is a simple text file containing
the instructions for how to generate each layer. Docker images are typically
quite large, often several GBs, while Dockerfiles are small and serve as
blueprints for the images. It is therefore good practice to have your
Dockerfile in your project Git repository, since it allows other users to
exactly replicate your project environment.

If you've been doing these tutorials on Windows you've been using the Docker
image `nbisweden/workshop-reproducible-research:slim`. The Dockerfile for
generating that image is called `Dockerfile_slim` and is located in your
`docker` directory (where you should hopefully be standing already). We will
now go through that file and discuss the different steps and what they do.
After that we'll build the image and test it out. Lastly, we'll start from that
image and make a new one to reproduce the results from the [Conda
tutorial](conda.md).

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
list of packages and reruns the second `RUN` command. **However, there is no
way for Docker to know that it should also update the apt-get package lists**.
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
with `docker images` that you can see your new image.

Now it's time to make our own Dockerfile to reproduce the results from the
[Conda tutorial](conda). If you haven't done the tutorial, it boils down to
creating a Conda environment file, setting up that environment, downloading
three RNA-seq data files, and running FastQC on those files. We will later
package and run the whole RNA-seq workflow in a Docker container, but for now
we keep it simple to reduce the size and time required.

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
can take a look below

??? example "Click to show an example of `Dockerfile_conda`"

    ```
    FROM my_docker_image:latest
    RUN conda config --add channels bioconda && \
        conda config --add channels conda-forge && \
        conda install -n base fastqc=0.11.9 sra-tools=2.10.1
    COPY run_qc.sh .
    CMD bash run_qc.sh
    ```

Build the image and tag it `my_docker_conda`:

```bash
docker build -t my_docker_conda -f Dockerfile_conda .
```

Verify that the image was built using `docker images`.

!!! note "Quick recap"
    In this section we've learned:

    * How the keywords `FROM`, `LABEL`, `MAINTAINER`, `RUN`, `ENV`, `SHELL`,
      `WORKDIR`, and `CMD` can be used when writing a Dockerfile.
    * The importance of letting each layer in the Dockerfile be a "logical
      unit".
    * How to use `docker build` to construct and tag an image from a Dockerfile.

## Managing containers

When you start a container with `docker run` it is given an unique id that you
can use for interacting with the container. Let's try to run a container from
the image we just created:

```bash
docker run my_docker_conda
```

If everything worked `run_qc.sh` is executed and will first download and then
analyse the three samples. Once it's finished you can list all containers,
including those that have exited.

```bash
docker container ls --all
```

This should show information about the container that we just ran. Similar to:
```
CONTAINER ID        IMAGE                    COMMAND                  CREATED             STATUS                      PORTS               NAMES
39548f30ce45        my_docker_conda     "/bin/bash -c 'bas..."    3 minutes ago       Exited (0) 3 minutes ago                             el
```

If we run `docker run` without any flags, your local terminal is attached to the
container. This enables you to see the output of `run_qc.sh`, but also disables
you from doing anything else in the meantime. We can start a container in
detached mode with the `-d` flag. Try this out and run `docker container ls`
to validate that the container is running.

By default, Docker keeps containers after they have exited. This can be
convenient for debugging or if you want to look at logs, but it also consumes
huge amounts of disk space. It's therefore a good idea to always run with
`--rm`, which will remove the container once it has exited.

If we want to enter a running container, there are two related commands we can
use, `docker attach` and `docker exec`. `docker attach` will attach local
standard input, output, and error streams to a running container. This can be
useful if your terminal closed down for some reason or if you started
a terminal in detached mode and changed your mind. `docker exec` can be used to
execute any command in a running container. It's typically used to peak in at
what is happening by opening up a new shell. Here we start the container in
detached mode and then start a new interactive shell so that we can see what
happens. If you use `ls` inside the container you can see how the script
generates file in the `data`, `intermediate` and `results` directories. Note
that you will be thrown out when the container exits, so you have to be quick.

```bash
docker run -d --rm --name my_container my_docker_conda
docker exec -it my_container /bin/bash
```

!!! tip
    Sometimes you would like to enter a stopped container. It's not a common
    use case, but it's included here for those of you who are doing these
    tutorials on Windows using Docker. Inadvertently shutting down your
    container can result in loss of a lot of work if you're not able to restart
    it. If you were to use `docker start` it would rerun the command set by
    `CMD`, which may not be what you want. Instead we use `docker commit
    container_name new_image_name` to convert the container `container_name` to
    the image `new_image_name`. We can then start a new container in that image
    as we normally would with `docker run -it --rm new_image_name /bin/bash`.
    Confusing, right? In theory, this would allow you to bypass using
    Dockerfiles and instead generate your image by entering an empty container
    in interactive mode, install everything there, and then commit as a new
    image. However, by doing this you would lose many of the advantages that
    Dockerfiles provide, such as easy distribution and efficient space usage
    via layers.

!!! note "Quick recap"
    In this section we've learned:

    * How to use `docker run` for starting a container and how the flags `-d`
      and `--rm` work.
    * How to use `docker container ls` for displaying information about the
      containers.
    * How to use `docker attach` and `docker exec` to interact with running
      containers.

### Bind mounts
There are obviously some advantages to isolating and running your data analysis
in containers, but at some point you need to be able to interact with the host
system to actually deliver the results. This is done via bind mounts. When you
use a bind mount, a file or directory on the *host machine* is mounted into
a container. That way, when the container generates a file in such a directory
it will appear in the mounted directory on your host system.

!!! tip
    Docker also has a more advanced way of data storage called
    [volumes](https://docs.docker.com/engine/admin/volumes/). Volumes provide
    added flexibility and are independent of the host machine's filesystem
    having a specific directory structure available. They are particularly
    useful when you want to share data *between* containers.

Say that we are interested in getting the resulting html reports from FastQC in
our container. We can do this by mounting a directory called, say,
`fastqc_results` in your current directory to the `/course/results/fastqc`
directory in the container. Try this out by running:

```bash
docker run --rm -v $(pwd)/fastqc_results:/course/results/fastqc my_docker_conda
```

Here the `-v` flag to docker run specifies the bind mount in the form of
`directory/on/your/computer:/directory/inside/container`. `$(pwd)` simply
evaluates to the working directory on your computer.

Once the container finishes validate that it worked by opening one of the html
reports under `fastqc_results/`.

We can also use bind mounts for getting files into the container rather than
out. We've mainly been discussing Docker in the context of packaging an
analysis pipeline to allow someone else to reproduce its outcome. Another
application is as a kind of very powerful environment manager, similarly to how
we've used Conda before. If you've organized your work into projects, then you
can mount the whole project directory in a container and use the container as
the terminal for running stuff while still using your normal OS for editing
files and so on. Let's try this out by mounting our current directory and start
an interactive terminal. Note that this will override the `CMD` command, so we
won't start the analysis automatically when we start the container.

```bash
docker run -it --rm -v $(pwd):/course/ my_docker_conda /bin/bash
```

If you run `ls` you will see that all the files in the `docker` directory are
there. Now edit `run_qc.sh` **on your host system** to download, say, 12000
reads instead of 15000. Then rerun the analysis with `bash run_qc.sh`. Tada!
Validate that the resulting html reports look fine and then exit the container
with `exit`.

## Distributing your images

There would be little point in going through all the trouble of making your
analyses reproducible if you can't distribute them to others. Luckily, sharing
Docker containers is extremely easy. The most common way is to use
[Dockerhub](https://hub.docker.com). Dockerhub lets you host unlimited public
images and one private image for free, after that they charge a small fee. If
you want to try it out here is how to do it:

1. Register for an account on [Dockerhub](https://hub.docker.com).

2. Use `docker login -u your_dockerhub_id` to login to the Dockerhub registry.

3. When you build an image, tag it with `-t your_dockerhub_id/image_name`,
   rather than just `image_name`.

4. Once the image has been built, upload it to Dockerhub with `docker push
   your_dockerhub_id/image_name`.

5. If another user runs `docker run your_dockerhub_id/image_name` the image
   will automatically be retrieved from Dockerhub. You can use `docker pull`
   for downloading without running.

If you want to refer to a Docker image in for example a publication, it's very
important that it's the correct version of the image. You can do this by adding
a tag to the name like this `docker build -t
your_dockerhub_id/image_name:tag_name`.

!!! tip
    On Dockerhub it is also possible to link to your Bitbucket or GitHub
    account and select repositories from which you want to automatically build
    and distribute Docker images. The Dockerhub servers will then build an
    image from the Dockerfile in your repository and make it available for
    download using `docker pull`. That way, you don't have to bother manually
    building and pushing using `docker push`. The GitHub repository for this
    course is linked to Dockerhub and the Docker images are built automatically
    from `Dockerfile` and `Dockerfile_slim`, triggered by changes made to the
    GitHub repository. You can take a look at the course on Dockerhub
    [here](https://hub.docker.com/r/nbisweden/workshop-reproducible-research).

## Packaging the case study

During these tutorials we have been working on a case study about the
multiresistant bacteria MRSA. Here we will build and run a Docker container
that contains all the work we've done so far.

* We've [set up a GitHub repository](git.md) for version control and for
  hosting our project.
* We've defined a [Conda environment](conda.md) that specifies the packages
  we're depending on in the project.
* We've constructed a [Snakemake workflow](snakemake.md) that performs the data
  analysis and keeps track of files and parameters.
* We've written a [R Markdown document](rmarkdown.md) that takes the results
  from the Snakemake workflow and summarizes them in a report.

The `docker` directory contains the final versions of all the files we've
generated in the other tutorials: `environment.yml`, `Snakefile`, `config.yml`,
`code/header.tex`, and `code/supplementary_material.Rmd`. The only difference
compared to the other tutorials is that we have also included the rendering of
the Supplementary Material HTML file into the Snakemake workflow as the rule
`make_supplementary`. Running all of these steps will take some time to execute
(around 20 minutes or so), in particular if you're on a slow internet
connection, and result in a 3.75 GB image.

Now take a look at `Dockerfile`. Everything should look quite familiar to you,
since it's basically the same steps as in the image we constructed in the
previous section, although some sections have been moved around. The main
difference is that we add the project files needed for executing the workflow
(mentioned in the previous paragraph), and install the conda packages listed in
`environment.yml`. If you look at the `CMD` command you can see that it will
run the whole Snakemake workflow by default.

Now run `docker build` as before, tag the image with `my_docker_project`:

````bash
docker build -t my_docker_project -f Dockerfile .
````
and go get a coffee while the image builds (or
you could use `docker pull nbisweden/workshop-reproducible-research` which
will download the same image).

Validate with `docker images`. Now all that remains is to run the whole thing
with `docker run`. We just want to get the results, so mount the directory
`/course/results/` to, say, `mrsa_results` in your current directory.

Well done! You now have an image that allows anyone to exactly reproduce your
analysis workflow (if you first `docker push` to Dockerhub that is).

!!! tip
    If you've done the [Jupyter Notebook tutorial](jupyter.md), you know that
    Jupyter Notebook runs as a web server. This makes it very well suited for
    running in a Docker container, since we can just expose the port Jupyter
    Notebook uses and redirect it to one of our own. You can then work with the
    notebooks in your browser just as you've done before, while it's actually
    running in the container. This means you could package your data, scripts
    and environment in a Docker image that also runs a Jupyter Notebook server.
    If you make this image available, say on Dockerhub, other researchers could
    then download it and interact with your data/code via the fancy interactive
    Jupyter notebooks that you have prepared for them. We haven't made any
    fancy notebooks for you, but we *have* set up a Jupyter Notebook server.
    Try it out if you want to (replace the image name with your version if
    you've built it yourself):

    ```bash
    docker run -it -p 8888:8888 nbisweden/workshop-reproducible-research \
        jupyter notebook  --ip=0.0.0.0 --allow-root
    ```

## Cleaning up

As mentioned before, Docker tends to consume a lot of disk space. In general,
`docker image rm` is used for removing images and `docker container rm` for
removing containers. Here are some convenient commands for cleaning up.

```bash
# Remove unused images
docker image prune

# Remove stopped containers
docker container prune

# Remove unused volumes (not used here, but included for reference)
docker volume prune

# Stop and remove ALL containers
docker container rm $(docker container ls -a -q)

# Remove ALL images
docker image rm $(docker image ls -a -q)
```
