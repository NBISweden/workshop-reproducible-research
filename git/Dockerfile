FROM rocker/verse:3.3.1
MAINTAINER "Rasmus Agren" rasmus.agren@scilifelab.se

WORKDIR /home
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
SHELL ["/bin/bash", "-c"]

# Install Miniconda3
RUN apt-get update && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends bzip2 curl
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -bf && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/root/miniconda3/bin:${PATH}"

#Install R packages
RUN install2.r --error \
    -r 'http://cran.rstudio.com' \
    -r 'http://bioconductor.org/packages/3.4/bioc' \
    ggplot2 \
    reshape2 \
    pheatmap \
    rtracklayer \
    GEOquery \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Set up the Conda environment
COPY environment.yml .
RUN conda env create -n my_conda_env -f environment.yml && \
    conda clean --all

# Add the workflow files
COPY Snakefile config.yml ./
COPY code ./code/

CMD source activate my_conda_env && snakemake --configfile config.yml
