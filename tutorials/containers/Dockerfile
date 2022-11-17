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
RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O && \
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
# Configure Conda channels and install Mamba
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install mamba \
    && mamba env update -n base -f environment.yml \
    && mamba install -c conda-forge jupyter \
    && conda clean --all

# Install sra-tools
RUN curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    && mkdir tmp_out && tar --wildcards -C tmp_out -zxvf sratoolkit.tar.gz */bin/* && cp -r tmp_out/*/bin/* $CONDA_PREFIX/bin/ \
    && rm -r sratoolkit.tar.gz tmp_out/

# Open port for running Jupyter Notebook
EXPOSE 8888

CMD snakemake -rp -c 1 --configfile config.yml
