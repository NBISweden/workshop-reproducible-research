FROM ubuntu:16.04

LABEL description = "Image for the NBIS reproducible research course."
MAINTAINER "John Sundh" john.sundh@scilifelab.se

# Use bash as shell
SHELL ["/bin/bash", "-c"]
# Set workdir
WORKDIR /course

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
RUN curl https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O && \
    bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/ && \
    rm Miniconda3-4.7.12.1-Linux-x86_64.sh && \
    /usr/miniconda3/bin/conda clean -tipsy && \
    ln -s /usr/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /usr/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Add conda to PATH and set locale
ENV PATH="/usr/miniconda3/bin:${PATH}"
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8

# Add project files
COPY environment.yml Snakefile config.yml ./
COPY code ./code/

# Install conda environment
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda env update -n base -f environment.yml && \
    conda clean --all

# Install jupyter and nb_conda
RUN conda install -c conda-forge jupyter nb_conda && \
    conda clean --all

# Install TinyTeX
ENV PATH="/usr/.TinyTeX/bin/x86_64-linux:${PATH}"
RUN echo -e "library(tinytex)\ntinytex::install_tinytex(dir='/usr/.TinyTeX')"| R --vanilla - && \
    tlmgr path add && \
    tlmgr update --self && \
    tlmgr install float grffile xcolor mdwtools epstopdf-pkg

# Open port for running Jupyter Notebook
# (Jupyter Notebook has to be separately installed in the container)
EXPOSE 8888

CMD snakemake -rp -c 1 --configfile config.yml
