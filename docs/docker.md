# Introduction to Docker
## What is Docker?
Docker is a tool designed to make it easier to create, deploy, and run applications by isolating them in "containers". The idea is to package your program together with everything it needs (other packages it depends on, various environment settings, data..) into one unit, i.e. a "container". This way we can ensure that the code generates exactly the same results regardless of where it's executed. Containers are in many ways similar to virtual machines, but more lightweight. Rather than creating a whole new OS they can use the same Linux kernel as the system that they're running on. While this might sound a bit technical, actually using Docker is quite easy, fun and very powerful.

Just as with Git, Docker was designed for software development but is rapidly becoming used also in scientific research. If you're doing web development you would for example run the webserver in one container and the database in another, thereby reducing the risk of one system affecting the other in unpredictable ways. Docker containers have also proven to be a very good solution to packaging, running and distributing scientific data analyses. Some applications relevant for reproducible research can be:

* When publishing, package your whole analysis pipeline together with your data in a Docker container and let it accompany the article. This way anyone can reproduce your analysis at the push of a button.
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

Here we show how to do it for Ubuntu, which is the most common desktop distribution. The same instructions apply to distributions based on Ubuntu, such as Elementary OS or Linux Mint. Docker requires a 64-bit Ubuntu version 14.04 or higher. If your OS is from 2015 or earlier you can double check this with `lsb_release -a`, if it's newer you're probably fine.

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

```
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

```
$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
ubuntu              latest              20c44cd7596f        2 weeks ago         123MB
```

We can now start a container running our image. We can refer to the image either by "REPOSITORY" or "IMAGE ID". The syntax for `docker run` is `docker run [OPTIONS] IMAGE [COMMAND] [ARG...]`. Let's run the command `uname -a` to get some info about the operating system. First run on your own system (skip this if you're using Windows, or use `ver` which is the Windows equivalent).

```
$ uname -a
Darwin liv433l.lan 15.6.0 Darwin Kernel Version 15.6.0: Mon Oct  2 22:20:08 PDT 2017; root:xnu-3248.71.4~1/RELEASE_X86_64 x86_64
```

Seems like I'm running the Darwin version of macOS. Then run it in the Ubuntu Docker image.

```bash
docker run ubuntu uname -a
```

And now I'm running on Linux! Try the same thing with `whoami`.

So, seems we can execute arbitary commands on Linux. Seems useful, but maybe a bit limited. We can also get an interactive terminal with the flags `-it`.

```
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

```
FROM ubuntu:16.04

LABEL description = "Minimal image for the NBIS reproducible research course."
MAINTAINER "Rasmus Agren" rasmus.agren@scilifelab.se
```

Here we use the instructions `FROM`, `LABEL` and `MAINTAINER`. The important one is `FROM`, which specifies which layer our image should start from. In this case we want it to be Ubuntu 16.04, which is one of the [official repositories](https://hub.docker.com/explore/). There are many roads to Rome when it comes to choosing the best image to start from. Say you want to run RStudio in a Conda environment through a Jupyter notebook. You could then start from one of the [rocker images](https://github.com/rocker-org/rocker) for R, a [Miniconda image](https://hub.docker.com/r/continuumio/miniconda/), or a [Jupyter image](https://hub.docker.com/r/jupyter/). Or you just start from one of the low-level official images and set up everything from scratch. `LABEL` and `MAINTAINER` is just meta-data that can be used for organising your various Docker components.

```
# Install Miniconda3 prerequisites, English language pack and fonts. Also include
# vim for convenience
RUN apt-get update && \
    apt-get install -y --no-install-recommends bzip2 wget language-pack-en fontconfig vim
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -bf && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/root/miniconda3/bin:${PATH}"
```

The next few lines introduce the important `RUN` instruction, which is used for executing shell commands. As a general rule you want each layer in an image to be a "logical unit". For example, if you want to install a program the `RUN` command should both retrive the program, install it and perform any necessary clean up. This is due to how layers work and how Docker decides what need's to be rerun between builds. The first command uses Ubuntu's package manager apt to install some packages (similar to how we've used Conda). Say that the first command was split into two instead:

```
# Update apt-get
RUN apt-get update

#Install packages
RUN apt-get install -y --no-install-recommends bzip2 wget language-pack-en fontconfig vim
```




Start from base depo rocker/biocontainers/jupyter

## Running containers

### Mounting volumes

## Distributing your Docker image

## Running the whole thing
