# HÆMA figure-generation runtime: publication-quality manuscript figures from endpoint tables.
# Base pinned by digest for reproducible builds (python:3.11-slim-bookworm, linux multi-arch) —
# the same digest used by haema-python so the two images share a verified base layer.
FROM python:3.11-slim-bookworm@sha256:e2d3af735aff6eeee600b1933bedd99da6645fedf572cc12ef4cc1331f2ceebe

LABEL org.opencontainers.image.title="haema-figures"
LABEL org.opencontainers.image.description="Matplotlib/seaborn + geopandas runtime for HÆMA publication figures"
LABEL org.opencontainers.image.version="0.4.0"

# matplotlib needs a writable config/cache dir when the image runs as a non-root UID
# (Nextflow uses -u uid:gid); without it the first import warns and can fail to cache fonts.
ENV MPLCONFIGDIR=/tmp/mplconfig \
    HOME=/tmp
RUN mkdir -p /tmp/mplconfig && chmod 1777 /tmp/mplconfig

# procps provides `ps`, which Nextflow requires inside the container to collect task metrics
# (the -slim base omits it).
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# Pinned visualisation stack. pandas/numpy are held at the same versions as haema-python so a
# table written by one image reads identically in the other.
RUN pip install --no-cache-dir \
    pandas==2.2.3 \
    numpy==2.1.3 \
    matplotlib==3.9.2 \
    seaborn==0.13.2

# Geospatial stack for the Ghana bioclimatic-zone map (Figure 1A). pyogrio/shapely/pyproj ship
# manylinux wheels that bundle GDAL/GEOS/PROJ, so no system GIS libraries are required.
RUN pip install --no-cache-dir \
    shapely==2.0.6 \
    pyproj==3.7.0 \
    pyogrio==0.10.0 \
    geopandas==1.0.1

# Fail the build early if the visualisation + geo stack cannot be imported.
RUN python -c "import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot, seaborn, pandas, numpy, geopandas, shapely, pyproj; print('haema-figures 0.4.0 deps OK')"
