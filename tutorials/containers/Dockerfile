FROM condaforge/miniforge3:24.7.1-0

LABEL authors="John Sundh, john.sundh@scilifelab.se; Erik Fasterius, erik.fasterius@nbis.se"
LABEL description="Image for the NBIS reproducible research course."

# Set workdir
WORKDIR /course

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Install required packages
RUN apt-get update && \
    apt-get install -y curl && \
    apt-get clean

# Install Quarto
ARG QUARTO_VERSION="1.3.450"
RUN mkdir -p /opt/quarto/${QUARTO_VERSION} && \
    curl -o quarto.tar.gz -L "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.tar.gz" && \
    tar -zxvf quarto.tar.gz -C "/opt/quarto/${QUARTO_VERSION}" --strip-components=1 && \
    rm quarto.tar.gz
ENV PATH=/opt/quarto/${QUARTO_VERSION}/bin:${PATH}

# Configure Conda
RUN conda config --set channel_priority strict && \
    conda config --append channels bioconda

# Install environment
COPY environment.yml ./
RUN conda env create -f environment.yml -n project_mrsa && \
    conda clean -a

RUN echo "source activate project_mrsa" >> ~/.bashrc
ENV PATH=/opt/conda/envs/project_mrsa/bin:${PATH}

# Add project files
COPY Snakefile config.yml ./
COPY code ./code/

CMD snakemake --configfile config.yml -p -c 1
