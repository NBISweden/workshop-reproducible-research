Containers can be large and complicated, but once you start using them regularly
you'll find that you start understand these complexities. There are lots of
different things you can do with images and containers in general, especially
when it comes to optimising build time or final image size. Here is some small
tips and tricks that you can be inspired from!

If you want to read more about containers in general you can check out these
resources:

* A "Get started with Docker" at [docker.com](https://docs.docker.com/get-started/).
* An [early paper](https://arxiv.org/abs/1410.0846) on the subject of using
  Docker for reproducible research.

## A base image with Conda

We've used Conda throughout this container tutorial, and we did it by
installing Conda inside the image when we built it. Wouldn't it be nice if we
didn't have to do this particular step? After all, installing Conda is just
busy-work, compared to installing the actual environment that we want to use
for the analyses. Luckily, there are already container images out there that
have Conda (and [Mamba](conda-4-extra-material)) installed, such as the ones
over at `condaforge/mambaforge`! What follows is a Dockerfile that you could
use instead of the ones described above to install things using a Conda
`environment.yml` file, without having to install Conda in the Docker image
when building it!

```no-highlight
FROM condaforge/mambaforge:4.10.1-0
LABEL description = "Image description"
MAINTAINER "Firstname Lastname" firstname.lastname@gmail.se

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /project

# Copy and install the Conda environment
COPY environment.yml ./
RUN conda config --set channel_priority strict \
    && mamba env update --name base --file environment.yml \
    && mamba clean --all --force-pkgs-dirs --yes

# Start Bash shell by default
CMD /bin/bash
```
