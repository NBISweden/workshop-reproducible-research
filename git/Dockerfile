FROM rocker/verse:3.4.0
MAINTAINER "Rasmus Agren" rasmus.agren@scilifelab.se

WORKDIR /home
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
SHELL ["/bin/bash", "-c"]

# Install Miniconda3
RUN apt-get update && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends bzip2 curl
RUN curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh --insecure -O && \
    bash Miniconda3-latest-Linux-x86_64.sh -bf && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/root/miniconda3/bin:${PATH}"

#Install R packages
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); \
    biocLite(suppressUpdates = T, pkgs = c('ggplot2', 'reshape2', 'pheatmap', 'rtracklayer', 'GEOquery'))" \
    && rm -rf /tmp/Rtmp* /tmp/*.rds

# Set up the Conda environment
COPY environment.yml .
RUN conda env update -n root -f environment.yml && \
    conda clean --all

# Add the workflow files
COPY Snakefile config.yml ./
COPY code ./code/

# Install jupyter and nb_conda
RUN conda install jupyter nb_conda

# Open port for running Jupyter Notebook
EXPOSE 8888

# Add a Jupyter Notebook profile
RUN mkdir -p -m 700 /root/.jupyter/ && \
    echo "c.NotebookApp.ip = '*'" >> /root/.jupyter/jupyter_notebook_config.py

CMD snakemake --configfile config.yml
