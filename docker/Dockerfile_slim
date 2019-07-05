FROM ubuntu:16.04

LABEL description = "Minimal image for the NBIS reproducible research course."
MAINTAINER "Leif Wigge" leif.wigge@scilifelab.se

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
                                               wget \
    && apt-get clean

# Install Miniconda and add to PATH
RUN curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O && \
    bash Miniconda3-4.6.14-Linux-x86_64.sh -bf -p /opt/miniconda3/ && \
    rm Miniconda3-4.6.14-Linux-x86_64.sh

# Add conda to PATH and set locale
ENV PATH="/opt/miniconda3/bin:${PATH}"
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /course

# Open port for running Jupyter Notebook
EXPOSE 8888

# Start Bash shell by default
CMD /bin/bash
