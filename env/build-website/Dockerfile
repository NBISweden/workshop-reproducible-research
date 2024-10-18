#
# First stage: Quarto installation
#
FROM ubuntu:20.04 AS quarto
ARG QUARTO_VERSION=1.4.549
ARG TARGETARCH
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
    && apt-get clean

RUN mkdir -p /opt/quarto \
    && curl -o quarto.tar.gz -L "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.tar.gz" \
    && tar -zxvf quarto.tar.gz -C /opt/quarto/ --strip-components=1 \
    && rm quarto.tar.gz

#
# Second stage: Conda environment
#
FROM condaforge/miniforge3:24.7.1-0
COPY --from=quarto /opt/quarto /opt/quarto
ENV PATH="${PATH}:/opt/quarto/bin"
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        pkg-config \
    && apt-get clean

# Copy and install the Conda environment
COPY environment.yml ./
RUN conda config --set channel_priority strict \
    && conda env update --name base --file environment.yml \
    && conda clean --all --force-pkgs-dirs --yes

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /work

# Start Bash shell by default
CMD /bin/bash

LABEL org.opencontainers.image.authors="erik.fasterius@nbis.se"
LABEL org.opencontainers.image.description="Dockerfile for the Tools for Reproducible Research NBIS / Elixir workshop"
