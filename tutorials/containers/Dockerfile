FROM condaforge/miniforge3:24.7.1-0

LABEL authors="John Sundh, john.sundh@scilifelab.se; Erik Fasterius, erik.fasterius@nbis.se"
LABEL description="Image for the NBIS reproducible research course."

# Set workdir
WORKDIR /course

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
ENV PATH /opt/quarto/${QUARTO_VERSION}/bin:${PATH}

# Configure Conda
RUN conda config --set channel_priority strict && \
    conda config --append channels bioconda

# Install environment
COPY environment.yml ./
RUN conda env update -f environment.yml -n base && \
    conda clean -a

# Add project files
COPY Snakefile config.yml ./
COPY code ./code/

# Open port for running Jupyter Notebook
EXPOSE 8888

CMD snakemake -p -c 1 --configfile config.yml
