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

# Configure Conda
RUN conda config --set channel_priority strict && \
    conda config --append channels bioconda

# Install environment
COPY environment.yml ./
RUN conda env create -f environment.yml -n project_mrsa && \
    conda clean -a

RUN echo "source activate project_mrsa" >> ~/.bashrc
ENV PATH=/opt/conda/envs/project_mrsa/bin:${PATH}
