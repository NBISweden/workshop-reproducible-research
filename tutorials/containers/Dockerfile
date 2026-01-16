FROM condaforge/miniforge3

LABEL authors="John Sundh, john.sundh@scilifelab.se; Erik Fasterius, erik.fasterius@nbis.se"
LABEL description="Minimal image for the NBIS reproducible research course."

WORKDIR /course

SHELL ["/bin/bash", "-c"]

# Install `curl` for downloading of FASTQ data later in the tutorial
RUN apt-get update && \
    apt-get install -y curl && \
    apt-get clean

# Configure Conda
RUN conda config --set channel_priority strict && \
    conda config --append channels bioconda

# Start Bash by default
CMD /bin/bash
