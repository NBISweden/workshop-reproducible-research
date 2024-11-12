FROM condaforge/miniforge3

LABEL authors="John Sundh, john.sundh@scilifelab.se"
LABEL description="Minimal image for the NBIS reproducible research course."

# Use bash as shell
SHELL ["/bin/bash", "--login", "-c"]

# Set workdir
WORKDIR /course

# Set timezone
ENV TZ="Europe/Stockholm"
ENV DEBIAN_FRONTEND=noninteractive

# Install package for setting timezone
RUN apt-get update && apt-get install -y tzdata curl && apt-get clean

# Configure Conda
RUN conda init bash && conda config --set channel_priority strict && \
    conda config --append channels bioconda && \
    conda config --append channels r && \
    conda config --set subdir linux-64

# Open port for running Jupyter Notebook
EXPOSE 8888

# Start Bash shell by default
CMD /bin/bash
