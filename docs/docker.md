# Introduction to Docker
## What is Docker?
Docker is a tool designed to make it easier to create, deploy, and run applications by isolating them in "containers". The idea is to package your program together with everything it needs (other packages it depends on, various environment settings, data..) into one unit, i.e. a "container". This way we can ensure that the code generates exactly the same results regardless of where it's executed. Containers are in many ways similar to virtual machines, but more lightweight. Rather than creating a whole new OS they can use the same Linux kernel as the system that they're running on. While this might sound a bit technical, actually using Docker is quite easy, fun and very powerful.

Just as with Git, Docker was designed for software development but is rapidly becoming used also in scientific research. If you're doing web development you would for example run the webserver in one container and the database in another, thereby reducing the risk of one system affecting the other in unpredictable ways. Docker containers have also proven to be a very good solution to packaging, running and distributing scientific data analyses. Some applications relevant for reproducible research can be:

* When publishing, package your whole analysis pipeline together with your data in a Docker image and let it accompany the article. This way anyone can reproduce your analysis at the push of a button.
* Packaging your analysis in a Docker container enables you to develop on e.g. your laptop and then seamlessly move to cluster or cloud to run the actual analysis.
* Say that you are collaborating on a project and you are using Mac and they are using Windows. Then you can set up a Docker container specific for your project to ensure that you're working in an identical environment.

All of this might sound a bit abstract so far, so let's get going.

## Tell me more
* A "Get started with Docker" at [docker.com](https://docs.docker.com/get-started/).

# Set up
This exercise depends on files from the course BitBucket repo. Take a look at the [intro](index) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/git_jupyter_docker`.

## Install Docker
First we need to install Docker. This is quite straightforward on macOS or Windows and a little more cumbersome on Linux. Note that Docker runs as root, which means that you have to have sudo privileges on your computer in order to install or run Docker.

### macOS
Go to  [https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac](https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac) and select "Get Docker for Mac (Stable)". This will download a dmg file. Click on it once it's done to start the installation. This will open up a window where you can drag the Docker.app to Applications. Close the window and click the Docker app from the Applications menu. Now it's basically just to click "next" a couple of times and we should be good to go. You can find the Docker icon in the menu bar in the upper right part of the screen.

### Windows
The instructions are different depending on if you have Windows 10 or Windows 7 (earlier versions aren't supported). In order to run Docker on Windows your computer must support Hardware Virtualization Technology and virtualization must be enabled. This is typically done in BIOS. This is outside the scope of this tutorial, so we'll simply go ahead as if though it's enabled and hope that it works.

On Windows 10 we will install Docker for Windows, which is available at [https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows](https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows). Select "Get Docker for Windows (Stable)".

1. Once it's downloaded, double-click Docker for Windows Installer.exe to run the installer.

2. Follow the install wizard and accept the license, authorize the installer, and proceed with the install. You will be asked to authorize Docker.app with your system password during the install process. Click Finish to exit the installer.

3. Start Docker from the Start menu. You can search for it if you cannot find it. The Docker whale icon should appear in the task bar.

4. Now we want to share your local drive(s), so that they are available for Docker. Right-click on the Docker whale icon in the task bar and select "Settings". Go to "Shared drives" and enable the drives you want Docker to have access to. Note that the drive where you'll be running the tutorials from has to be enabled (most likely `C:\`).

On Windows 7 we will instead use Docker Toolbox, which is available at [https://docs.docker.com/toolbox/toolbox_install_windows/](https://docs.docker.com/toolbox/toolbox_install_windows/). Select "Get Docker Toolbox for Windows".

1. Install Docker Toolbox by double-clicking the installer. Step through the installation and accept all the defaults. If Windows security dialog prompts you to allow the program to make a change, choose Yes. If you get a promt asking "Would you like to install this device software?" select Install.

2. You should now have a Docker Quickstart icon on the desktop.

TODO: FIX VOLUME SHARING. TALK ABOUT /C/

### Linux
How to install Docker differs a bit depending on your Linux distro, but the steps are the same. For details on how to do it on your distro see [https://docs.docker.com/engine/installation/#server](https://docs.docker.com/engine/installation/#server).

Here we show how to do it for Ubuntu, which is the most common desktop distribution. The same instructions apply to distributions based on Ubuntu, such as Elementary OS or Linux Mint. Docker requires a 64-bit Ubuntu version 14.04 or higher. If your OS is from 2015 or earlier you can double check this with `lsb_release -a`. If it's newer you're probably fine.

1. Add the GPG key for the official Docker repository to the system:
```bash
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
```

2. Add the Docker repository to APT sources:
```bash
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
```

3. Update the package database with the Docker packages from the new repo:
```bash
sudo apt-get update
```

4. Install Docker Community Edition:
```bash
sudo apt-get install -y docker-ce
```

5. Docker should now be installed, the daemon started, and the process enabled to start on boot. Check that it's running:
```bash
sudo systemctl status docker
```
The output should say something about "Active: active (running) since..".

As mentioned before, Docker needs to run as root. You can archive this by prepending all Docker commands with `sudo`. This is the approch that we will take in this tutorial, since the set up becomes a little simpler. If you plan on continuing using Docker you can get rid of this by adding your user to the group `docker`. Here are instructions for how to do this: [https://docs.docker.com/engine/installation/linux/linux-postinstall/](https://docs.docker.com/engine/installation/linux/linux-postinstall/).

# Practical exercise
## The very basics
We're almost ready to start, just one last note on nomenclature. You might have noticed that we sometimes refer to "Docker images" and sometimes to "Docker containers". A container is simply an instance of an image. You can have an image containing, say, a certain Linux distribution, and then start multiple containers running that same OS.

> NOTE: If you don't have root privileges you have to prepend all Docker commands with `sudo ` (see Linux set up instructions)

Docker containers typically run Linux, so let's start by downloading an image containing Ubuntu (a popular Linux distribution that is based on only open-source tools).

```no-highlight
$ docker pull ubuntu:latest
Using default tag: latest
latest: Pulling from library/ubuntu
660c48dd555d: Pull complete
4c7380416e78: Pull complete
421e436b5f80: Pull complete
e4ce6c3651b3: Pull complete
be588e74bd34: Pull complete
Digest: sha256:7c67a2206d3c04703e5c23518707bdd4916c057562dd51c74b99b2ba26af0f79
Status: Downloaded newer image for ubuntu:latest
```

You might have noticed that it downloaded five different layers with weird hashes as names ("660c48dd555d" and so on). This represents a very fundamental property of Docker images that we'll get back to in just a little while. For now let's just look at our new collection of Docker images:

```no-highlight
$ docker image ls
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
ubuntu              latest              20c44cd7596f        2 weeks ago         123MB
```

We can now start a container running our image. We can refer to the image either by "REPOSITORY" or "IMAGE ID". The syntax for `docker run` is `docker run [OPTIONS] IMAGE [COMMAND] [ARG...]`. Let's run the command `uname -a` to get some info about the operating system. First run on your own system (skip this if you're using Windows, or use `ver` which is the Windows equivalent).

```no-highlight
$ uname -a
Darwin liv433l.lan 15.6.0 Darwin Kernel Version 15.6.0: Mon Oct  2 22:20:08 PDT 2017; root:xnu-3248.71.4~1/RELEASE_X86_64 x86_64
```

Seems like I'm running the Darwin version of macOS. Then run it in the Ubuntu Docker image.

```bash
docker run ubuntu uname -a
```

And now I'm running on Linux! Try the same thing with `whoami`.

So, seems we can execute arbitary commands on Linux. Seems useful, but maybe a bit limited. We can also get an interactive terminal with the flags `-it`.

```no-highlight
$ docker run -it ubuntu
root@1f339e929fa9:/#
```

Here you can do whatever; install, run, remove stuff; it will still be within the container and never affect your host system. Now exit the container with `exit`.

Ok, so Docker let's us work in any OS in a quite convenient way. That would probably be useful on its own, but it's much more powerful than that.

## Building a Docker image
In the previous section we downloaded a Docker image of Ubuntu and noticed that it was based on layers, each with a unique hash as id. An image in Docker is based on a number of read-only layers, where each layer contains the differences to the previous layers. If you've done the [Git tutorial](git) this might remind you of how a Git commit contains the difference to the previous commit. The great thing about this is that we can start from one base layer, say containing an operating system and some utility programs, and then generate many new images based on this, say 10 different project-specific images. The total space requirements would then only be `base_image + 10*(project_specific)` rather than `10*(base_image + project_specific)`. For example, Bioconda (see the [Conda tutorial](conda)) has one base image and then one layer for each of the >3000 packages available in Bioconda.

Docker provides a convenient way to describe how to go from a base image to the image we want by using a "Dockerfile". This is a simple text file containing the instructions for how to generate each layer. Docker images are typically quite large, often several GBs, while Dockerfiles are small and serve as blueprints for the images. It is therefore good practice to have a Dockerfile in your project Git repository, since it allows other users to exactly replicate your project environment.

If you've been doing these tutorials on Windows you've been using the Docker image `scilifelablts/reproducible_research_course_slim`. The Dockerfile for generating that image is called `Dockerfile_slim` and is located in your `git_jupyter_docker` directory. We will now go though that file and discuss the different steps and what they do. After that we'll build the image and test it out. Lastly, we'll start from that image and make a new one to reproduce the results from the [Conda tutorial](conda).

Here are the first few lines of `Dockerfile_slim`. Each line in the Dockerfile will typically result in one layer in the resulting image. The format for Dockerfiles is `INSTRUCTION arguments`. A full specification of the format, together with best practices, can be found [here](https://docs.docker.com/engine/reference/builder/).

```no-highlight
FROM ubuntu:16.04

LABEL description = "Minimal image for the NBIS reproducible research course."
MAINTAINER "Rasmus Agren" rasmus.agren@scilifelab.se
```

Here we use the instructions `FROM`, `LABEL` and `MAINTAINER`. The important one is `FROM`, which specifies which layer our image should start from. In this case we want it to be Ubuntu 16.04, which is one of the [official repositories](https://hub.docker.com/explore/). There are many roads to Rome when it comes to choosing the best image to start from. Say you want to run RStudio in a Conda environment through a Jupyter notebook. You could then start from one of the [rocker images](https://github.com/rocker-org/rocker) for R, a [Miniconda image](https://hub.docker.com/r/continuumio/miniconda/), or a [Jupyter image](https://hub.docker.com/r/jupyter/). Or you just start from one of the low-level official images and set up everything from scratch. `LABEL` and `MAINTAINER` is just meta-data that can be used for organising your various Docker components.

```no-highlight
# Install Miniconda3 prerequisites, English language pack and fonts. Also include
# vim for convenience
RUN apt-get update && \
    apt-get install -y --no-install-recommends bzip2 wget language-pack-en fontconfig vim
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -bf && \
    rm Miniconda3-latest-Linux-x86_64.sh
```

The next few lines introduce the important `RUN` instruction, which is used for executing shell commands. As a general rule you want each layer in an image to be a "logical unit". For example, if you want to install a program the `RUN` command should both retrive the program, install it and perform any necessary clean up. This is due to how layers work and how Docker decides what need's to be rerun between builds. The first command uses Ubuntu's package manager apt to install some packages (similar to how we've previously used Conda). Say that the first command was split into two instead:

```no-highlight
# Update apt-get
RUN apt-get update

#Install packages
RUN apt-get install -y --no-install-recommends bzip2 wget language-pack-en fontconfig vim
```

The first command will update the apt-get package lists and the seconds will install the packages bzip2, wget, language-pack-en, fontconfig, and vim. Say that you build this image now, and then in a month's time you realize that you would have liked a Swedish language pack instead of an English. You change to `language-pack-sv` and rebuild the image. Docker detects that there is no layer with the new list of packages and re-runs the second `RUN` command. *However, there is no way for Docker to know that it should also update the apt-get package lists*. You therefore risk to end up with old versions of packages and, even worse, the versions would depend on when the previous version of the image was first built.

The next `RUN` command retrieves and installs Miniconda3. Let's see what would happen if we had that as three separate commands instead.

```no-highlight
# Download Miniconda3
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install it
RUN bash Miniconda3-latest-Linux-x86_64.sh -bf

# Remove the downloaded installation file
RUN rm Miniconda3-latest-Linux-x86_64.sh
```

Remember that each layer contains the difference compared to the previous layer? What will happen here is that the first command adds the installation file and the second will unpack the file and install the software. The third layer will say "the installation file should no longer exist". However, the file will still remain in the image since the image is constructed layer-by-layer bottom-up. This results in unnecessarily many layers and bloated images.

```no-highlight
# Add conda to PATH and set locale
ENV PATH="/root/miniconda3/bin:${PATH}"
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
```

Here we use the new instruction `ENV`. The first command adds `conda` to the path, so we can write `conda install` instead of `/root/miniconda3/bin/conda install`. The other two set an UTF-8  character encoding so that we can use weird characters (and a bunch of other things).

```no-highlight
# Install git, nano and ca-certificates from conda-forge
RUN conda install -c conda-forge git ca-certificates nano && \
    conda clean --all
```

Here we install some packages with Conda. Note that we run `conda clean --all` to remove any downloaded packages afterwards.

```no-highlight
# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /home

# Start Bash shell by default
CMD /bin/bash
```

`SHELL` simply sets which shell to use and `WORKDIR` which directory the container should start in. `CMD` is an interesting instruction. It sets what a container should run when nothing else is specified. It can be used for example for printing some information on how to use the image or, as here, start a shell for the user. If the purpose of your image is to accompany a publication then `CMD` could be to run the workflow that generates the paper figures from raw data.

Ok, so now we understand how a Dockerfile works. Constructing the image from the Dockerfile is really simple. Try it out now.

```no-highlight
$ docker build -f Dockerfile_slim -t my_docker_image .
Step 1/12 : FROM ubuntu:16.04
 ---> 20c44cd7596f
Step 2/12 : LABEL description = "Minimal image for the NBIS reproducible research course."
.
.
[lots of stuff]
.
.
Step 12/12 : CMD "/bin/bash"
 ---> Running in f34c2dbbbecf
 ---> aaa39bdeb78a
Removing intermediate container f34c2dbbbecf
Successfully built aaa39bdeb78a
Successfully tagged my_docker_image:latest
```

`-f` sets which Dockerfile to use and `-t` tags the image with a name. This name is how you will refer to the image later. Lastly, the `.` is the path where the image should be build (`.` means the current directory). This had no real impact in this case, but matters if you want to import files. Validate with `docker image ls` that you can see your new image.

Now it's time to make our own Dockerfile to reproduce the results from the [Conda tutorial](conda). If you haven't done the tutorial, it boils down to creating a Conda environment file, setting up that environment, downloading three RNA-seq data files, and running FastQC on those files. We will later package and run the whole RNA-seq workflow in a Docker container, but for now we keep it simple to reduce the size and time required.

The Conda tutorial uses a shell script, `run_qc.sh`, for downloading and running the analysis. If we want to use the same script we need to include it in the image. To do that we first copy it to our current directory.

```bash
cp ../conda/code/run_qc.sh .
```

So, this is what we need to do:

0. Create the file `Dockerfile_conda`.
1. Set `FROM` to the image we just built.
2. Install the required packages with Conda. We could do this by adding `environment.yml` from the Conda tutorial, but here we do it directly as `RUN` commands. We need the add the conda-forge and bioconda channels with `conda config --add channels channel_name` and install `fastqc=0.11` and `sra-tools=2.8` with `conda install`. There is little point in defining and activating a Conda environment since the container is self-contained, but do so if you want.
3. Add `run_qc.sh` to the image by using the `COPY` instruction. The syntax is `COPY source target`, so in our case simply `COPY run_qc.sh .` to copy to the home directory in the image.
4. Set the default command for the image to `bash run_qc.sh`, which will execute the shell script.
5. Build the image and tag it `my_docker_conda`. Verify with `docker image ls`.

## Managing containers
When you start a container with `docker run` it is given an unique id that you can use for interacting with the container. Let's try to run the image we just created:

```bash
docker run my_docker_conda
```

If everything worked `run_qc.sh` is executed and will first download and then analyse the three samples. Once it's finished you can list all containers, including those that have exited.

```no-highlight
$ docker container ls --all
CONTAINER ID        IMAGE               COMMAND                   CREATED             STATUS                           PORTS               NA
MES
39548f30ce45        my_docker_conda     "/bin/bash -c 'bas..."    3 minutes ago       Exited (0) 3 minutes ago                             el
```

If we run `docker run` without any flags, your local terminal is attached to the container. This enables you to see the output of `run_qc.sh`, but also disables you from doing anything else in the meantime. We can start a container in detached mode with the `-d` flag. Try this out and run `docker container ls` to validate that the container is running.

By default Docker keeps containers after they have exited. This can be convenient for debugging or if you want to look at logs, but it also consumes huge amounts of disk space. It's therefore a good idea to always run with `--rm`, which will remove the container once it has exited.

If we want to enter a running container, there are two related commands we can use, `docker attach` and `docker exec`. `docker attach` will attach local standard input, output, and error streams to a running container. This can be useful if your terminal closed down for some reason or if you started a terminal in detached mode and changed your mind. `docker exec` can be used to execute any command in a running container. It's typically used to peak in at what is happening by opening up a new shell. Here we start the container in detached mode and then start a new interactive shell so that we can see what happens. If you use `ls` inside the container you can see how the script generates file in the `data`, `intermediate` and `results` directories. Note that you will be thrown out when the container exits though.

```bash
docker run -d --rm --name my_container my_docker_conda
docker exec -it my_container /bin/bash
```

Sometimes you would like to enter a stopped container. It's not a common use case, but it's included here for those of you who are doing these tutorials on Windows using Docker. If you for some reason happen to shut down your container it can result in loss of a lot of work if you're not able to restart it. If you were to use `docker start` it would rerun the command set by `CMD`, which may not be what you want. Instead we use `docker commit container_name new_image_name` to convert the container `container_name` to the image `new_image_name`. We can the start a new container in that image as we normally would with `docker run -it --rm new_image_name`. Weird, right? In theory this would allow you to bypass using Dockerfiles and instead generate your image by entering an empty container in interactive mode, install everything there, and then commit as a new image. However, by doing this you would lose many of the advantages that Dockerfiles provide, such as easy distribution and efficient space usage via layers.

### Bind mounts
There are obviously some advantages to isolating and running a data analysis in containers, but at some point you need to be able to interact with the host system to actually deliver the results. This is done via bind mounts. When you use a bind mount, a file or directory on the *host machine* is mounted into a container. That way, when the container generates a file in such a directory it will appear in the mounted directory on your host system.

!!! tip
    Docker also has a more advanced way of data storage called [volumes](https://docs.docker.com/engine/admin/volumes/). Volumes provide added flexibility and are independent of the host machine's filesystem having a specific directory structure available. They are particularly useful when you want to share data *between* containers.

Say that we are interested in getting the resulting html reports from FastQC in our container. We can do this by mounting a directory called, say, `fastqc_results` in your current directory to the `/home/results/fastqc` directory in the container. Validate that it worked by opening one of the html reports.

```bash
docker run --rm -v $(pwd)/fastqc_results:/home/results/fastqc my_docker_conda
```

We can also use bind mounts for getting files into the container rather than out. We've mainly been discussing Docker in the context of  packaging an analysis pipeline to allow anyone to reproduce its outcome. Another application is as a kind of very powerful environment manager, similarly to how we've used Conda before. If you've organised your work into projects, then you can mount the whole project directory in a container and use the container as the terminal for running stuff while still using your normal OS for editing files and so on. Let's try this out by mounting our current directory and start an interactive terminal. Note that this will override the `CMD` command, so we won't start the analysis automatically when we start the container.

TODO: Validate this on Windows

```bash
docker run -it --rm -v $(pwd):/home/ my_docker_conda /bin/bash
```

If you run `ls` you will see that all the files in the `git_jupyter_docker` directory are there. Now edit `run_qc.sh` *on your host system* to download, say, 15000 reads instead of 12000. Then rerun the analysis with `bash run_qc.sh`. Tada! Validate that the resulting html reports look fine and then exit the container with `exit`.

## Distributing your images
There would be little point in going through all the trouble of making your analyses reproducible if you can't distribute them to others. Luckily, sharing Docker containers is extremely easy. The most common way is to use [Dockerhub](https://hub.docker.com). Dockerhub lets you host unlimited public images and one private image for free, after that they charge a small fee. If you want to try it out here is how to do it:

1. Register for an account on [Dockerhub](https://hub.docker.com).
2. Use `docker login -u your_dockerhub_id` to login to the Dockerhub registry.
3. When you build an image, tag it with `-t your_dockerhub_id/image_name`, rather than just `image_name`.
4. Once the image has been built, upload it to Dockerhub with `docker push your_dockerhub_id/image_name`.
5. If another user runs `docker run your_dockerhub_id/image_name` the image will automatically be retrieved from Dockerhub. You can use `docker pull` for downloading without running.

That was easy!

If you want to refer to a Docker image in for example a publication, it's very important that it's the correct version of the image. You can do this by adding a tag to the name like this `docker build -t your_dockerhub_id/image_name:tag_name`.

## Packaging and running the MRSA workflow
During these tutorials we have been working on a case study about the multiresistant bacteria MRSA. Here we will build and run a Docker container that contains all the work we've done. This will take some time to execute (~20 min or so), in particular if you're on a slow internet connection, and result in a 5.4 GB image.

* We've [set up a Bitbucket repository](git) for version control and for hosting our project.
* We've defined a [Conda environment](conda) that specifies the packages we're depending on in the project.
* We've constructed a [Snakemake workflow](snakemake) that performs the data analysis and keeps track of files and parameters.
* We've written a [R Markdown document](rmarkdown) that takes the results from the Snakemake workflow and summarizes them in a report.

The `git_jupyter_docker` directory contains the final versions of all the files we've generated in the other tutorials: `environment.yml`, `Snakefile`, `config.yml`, `code/header.tex`, and `code/supplementary_material.Rmd`. The only difference compared to the tutorials is that we have also included the rendering of the supplementary material pdf into the Snakemake workflow as the rule `make_supplementary`.

Now take a look at `Dockerfile`. Everything should look quite familiar to you, since it's basically the same steps as in the image we constructed in the previous section. The main difference is that here we start from `rocker/verse:3.3.1` rather than from `ubuntu:16.04`. This image contains RStudio and a number of publishing-related packages, most notably LaTeX for generating PDF reports. We need this in order to be able to render the Supplementary material report to PDF, but it also takes up quite a lot of space (2.17 GB). The other main difference is that we install a number of R packages with `install2.r`. We could have included these in the Conda environment as well, but R's package management is quite good so we might as well use that. If you look at the `CMD` command you can see that it will activate the Conda environment and run the whole Snakemake workflow by default.

Now run `docker build` as before and go get a coffe while the image builds. Validate with `docker image ls`. Now all that remains is to run the whole thing with `docker run`. We just want to get the results, so mount the directory `/home/results/` to, say, `mrsa_results` in your current directory. Well done! You now have an image that allows anyone to exactly reproduce your analysis workflow (if you first `docker push` to Dockerhub that is).

## Cleaning up
As mentioned before, Docker tends to consume a lot of disk space. In general, `docker image rm` is used for removing images and `docker container rm` for removing containers. Here are some convenient commands for cleaning up.

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
