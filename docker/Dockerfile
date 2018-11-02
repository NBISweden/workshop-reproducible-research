FROM rocker/verse:3.4.0
MAINTAINER "Rasmus Agren" rasmus.agren@scilifelab.se

WORKDIR /course
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
SHELL ["/bin/bash", "-c"]

# Install Miniconda3
RUN apt-get update && apt-get -y dist-upgrade
RUN apt-get install -y --no-install-recommends bzip2 curl ca-certificates
RUN curl https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O && \
    bash Miniconda3-4.5.11-Linux-x86_64.sh -bf -p /opt/miniconda3/ && \
    rm Miniconda3-4.5.11-Linux-x86_64.sh

# Add Conda to PATH
ENV PATH="/opt/miniconda3/bin:${PATH}"

# Install R packages
RUN Rscript -e "BiocInstaller::biocLite(suppressUpdates = T, \
    pkgs = c('ggplot2', 'reshape2', 'pheatmap', 'rtracklayer', 'devtools'))" \
    && rm -rf /tmp/Rtmp* /tmp/*.rds
# Install devel version due to https://github.com/seandavi/GEOquery/issues/39
RUN Rscript -e "devtools::install_github('seandavi/GEOquery')" \
    && rm -rf /tmp/Rtmp* /tmp/*.rds

# Set up the Conda environment
COPY environment.yml .
RUN conda env update -n base -f environment.yml && \
    conda clean --all

# Add the workflow files
COPY Snakefile config.yml ./
COPY code ./code/

# Install jupyter and nb_conda
RUN conda install -c conda-forge jupyter nb_conda && \
    conda clean --all

# Open port for running Jupyter Notebook
EXPOSE 8888

# Add a Jupyter Notebook profile
RUN mkdir -p -m 700 /root/.jupyter/ && \
    echo "c.NotebookApp.ip = '0.0.0.0'" >> /root/.jupyter/jupyter_notebook_config.py

CMD snakemake --configfile config.yml
