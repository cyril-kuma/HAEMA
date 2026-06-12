# HÆMA R/Bioconductor runtime for formal phyloseq + decontam endpoints.
# Base pinned by digest for reproducible builds (rocker/r-ver:4.4.2 -> Bioconductor 3.20).
FROM rocker/r-ver:4.4.2@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c

LABEL org.opencontainers.image.title="haema-r"
LABEL org.opencontainers.image.description="R/Bioconductor runtime for HÆMA phyloseq and decontam outputs"
LABEL org.opencontainers.image.version="0.3.0"

# procps provides `ps`, which Nextflow requires inside the container to collect task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
        libhdf5-dev \
        procps \
    && rm -rf /var/lib/apt/lists/*

# build_r_outputs.R only uses base R plus phyloseq:: and decontam::; ape is pulled in
# transitively by phyloseq and is verified by R_BIOCONDUCTOR_PREFLIGHT. tidyverse was
# previously installed but is unused, so it is dropped to cut build time and image size.
# Bioconductor release is pinned to 3.20 (the release mapped to R 4.4.x) for reproducibility.
RUN Rscript -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)' \
    && Rscript -e 'BiocManager::install(c("phyloseq", "decontam", "ape"), ask = FALSE, update = FALSE)'

# Fail the build early if the formal Bioconductor stack is not importable.
RUN Rscript -e 'stopifnot(all(vapply(c("phyloseq","decontam","ape"), requireNamespace, logical(1), quietly = TRUE))); cat("haema-r deps OK\n")'
