FROM ubuntu:16.04

LABEL description = "Image for the NBIS reproducible research course."
MAINTAINER "Leif Wigge" leif.wigge@scilifelab.se

WORKDIR /course
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
ENV HOME /usr
SHELL ["/bin/bash", "-c"]

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

# Install TinyTeX
RUN wget -qO- "https://yihui.name/gh/tinytex/tools/install-unx.sh" | sh  && \
    tlmgr update --self && \
    tlmgr install float && \
    tlmgr install grffile

# Install Miniconda and add to PATH
RUN curl https://repo.continuum.io/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O && \
    bash Miniconda3-4.6.14-Linux-x86_64.sh -bf -p /usr/miniconda3/ && \
    rm Miniconda3-4.6.14-Linux-x86_64.sh
ENV PATH="/usr/miniconda3/bin:${PATH}"

# Add project files
COPY environment.yml Snakefile config.yml ./
COPY code ./code/

# Install conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda env update -n base -f environment.yml && \
    conda clean --all

# Install jupyter and nb_conda
RUN conda install -c conda-forge jupyter nb_conda && \
    conda clean --all

# Open port for running Jupyter Notebook
# (Juoyter Notebook has to be separately installed in the container)
EXPOSE 8888

CMD snakemake -rp --configfile config.yml
