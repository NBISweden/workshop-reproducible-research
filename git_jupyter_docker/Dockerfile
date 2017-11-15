FROM rocker/verse:3.3.1
MAINTAINER "Rasmus Agren" rasmus.agren@scilifelab.se

RUN apt-get update && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends bzip2

WORKDIR /home
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
SHELL ["/bin/bash", "-c"]

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -bf && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/root/miniconda3/bin:${PATH}"

ADD environment.yml .
RUN conda env create -n my_conda_env -f environment.yml

ADD Snakefile .
CMD source activate my_conda_env && /bin/bash
