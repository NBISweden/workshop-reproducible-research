#
# First stage: Quarto installation
#
FROM ubuntu:20.04 AS quarto
ARG QUARTO_VERSION=1.10.18
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
ARG TARGETARCH
COPY --from=quarto /opt/quarto /opt/quarto
ENV PATH="${PATH}:/opt/quarto/bin"
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        pkg-config \
    && apt-get clean

# Install Chrome Headless Shell, which Quarto uses to render Mermaid diagrams
# to SVG at build time (only available for amd64)
RUN if [ "${TARGETARCH}" = "amd64" ]; then \
        apt-get update \
        && apt-get install -y --no-install-recommends \
            fonts-dejavu-core \
            libasound2 \
            libatk-bridge2.0-0 \
            libatk1.0-0 \
            libcups2 \
            libdbus-1-3 \
            libdrm2 \
            libgbm1 \
            libnss3 \
            libpango-1.0-0 \
            libxcomposite1 \
            libxdamage1 \
            libxfixes3 \
            libxkbcommon0 \
            libxrandr2 \
            unzip \
        && apt-get clean \
        && quarto install chrome-headless-shell --no-prompt; \
    fi

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
LABEL org.opencontainers.image.description="Dockerfile for building the Tools for Reproducible Research NBIS / Elixir workshop website"
