# HÆMA Python scientific runtime for UMAP/HDBSCAN mixed-template denoising.
# Base pinned by digest for reproducible builds (python:3.11-slim-bookworm, linux multi-arch).
FROM python:3.11-slim-bookworm@sha256:e2d3af735aff6eeee600b1933bedd99da6645fedf572cc12ef4cc1331f2ceebe

LABEL org.opencontainers.image.title="haema-python"
LABEL org.opencontainers.image.description="Python runtime for HÆMA ONT blood-meal metabarcoding (UMAP/HDBSCAN denoising)"
LABEL org.opencontainers.image.version="0.3.0"

# numba (used by UMAP) needs a writable cache + HOME when the image runs as a non-root UID
# (Nextflow uses -u uid:gid). Default to writable /tmp locations; pipeline processes override
# NUMBA_CACHE_DIR to the task work dir. Without this, umap import fails and denoising falls back.
ENV NUMBA_CACHE_DIR=/tmp/numba_cache \
    HOME=/tmp
RUN mkdir -p /tmp/numba_cache && chmod 1777 /tmp/numba_cache

# procps provides `ps`, which Nextflow requires inside the container to collect task metrics
# (the -slim base omits it). build-essential/gfortran/openblas are needed to build hdbscan.
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential gfortran libopenblas-dev liblapack-dev procps \
    && rm -rf /var/lib/apt/lists/*

# Pinned scientific stack. numpy/scipy/scikit-learn are transitive UMAP/HDBSCAN
# dependencies but are pinned explicitly so the resolved build is reproducible.
RUN pip install --no-cache-dir \
    biopython==1.84 \
    pandas==2.2.3 \
    numpy==2.1.3 \
    scipy==1.14.1 \
    scikit-learn==1.5.2 \
    umap-learn==0.5.7 \
    hdbscan==0.8.40

# Fail the build early if the scientific denoising stack cannot be imported.
RUN python -c "import numpy, scipy, sklearn, umap, hdbscan, Bio, pandas; print('haema-python deps OK')"
