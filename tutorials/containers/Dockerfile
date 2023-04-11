FROM condaforge/mambaforge

LABEL description = "Image for the NBIS reproducible research course."
MAINTAINER "John Sundh" john.sundh@scilifelab.se

# Use bash as shell
SHELL ["/bin/bash", "--login", "-c"]

# Set workdir
WORKDIR /course

# Configure Conda/Mamba
RUN mamba init bash && conda config --set channel_priority strict && \
    conda config --append channels bioconda && \
    conda config --append channels r && \
    conda config --set subdir linux-64

# Install environment
COPY environment.yml ./
RUN mamba env create -f environment.yml -n project_mrsa && \
    mamba clean -a

# Set mrsa-workflow environment as active at start-up
RUN echo "source activate project_mrsa" >> ~/.bashrc

# Add environment to PATH
ENV PATH /opt/conda/envs/project_mrsa/bin:${PATH}

# Add project files
COPY Snakefile config.yml ./
COPY code ./code/

# Open up port 8888
EXPOSE 8888

CMD snakemake -rp -c 1 --configfile config.yml