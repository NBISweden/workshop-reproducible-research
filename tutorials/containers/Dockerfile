FROM --platform=amd64 condaforge/miniforge3

LABEL description="Image for the NBIS reproducible research course."
LABEL author="John Sundh"
LABEL email="john.sundh@scilifelab.se"

# Use bash as shell
SHELL ["/bin/bash", "--login", "-c"]

# Set workdir
WORKDIR /course

# Set timezone
ENV TZ="Europe/Stockholm"
ENV DEBIAN_FRONTEND=noninteractive

# Install packages require for timezone and Quarto installation
RUN apt-get update \
    && apt-get install -y tzdata curl \
    && apt-get clean

# Install Quarto
ARG QUARTO_VERSION="1.3.450"
RUN mkdir -p /opt/quarto/${QUARTO_VERSION} \
    && curl -o quarto.tar.gz -L "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.tar.gz" \
    && tar -zxvf quarto.tar.gz -C "/opt/quarto/${QUARTO_VERSION}" --strip-components=1 \
    && rm quarto.tar.gz
ENV PATH /opt/quarto/${QUARTO_VERSION}/bin:${PATH}

# Configure Conda
RUN conda init bash && conda config --set channel_priority strict && \
    conda config --append channels bioconda && \
    conda config --append channels r && \
    conda config --set subdir linux-64

# Install environment
COPY environment.yml ./
RUN conda env create -f environment.yml -n project_mrsa && \
    conda clean -a

# Set mrsa-workflow environment as active at start-up
RUN echo "source activate project_mrsa" >> ~/.bashrc

# Add environment to PATH
ENV PATH /opt/conda/envs/project_mrsa/bin:${PATH}

# Add project files
COPY Snakefile config.yml ./
COPY code ./code/

# Open up port 8888
EXPOSE 8888

CMD snakemake -p -c 1 --configfile config.yml
