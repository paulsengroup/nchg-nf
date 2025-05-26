# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.13-slim AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&& apt-get install -y procps \
&& rm -rf /var/lib/apt/lists/*

RUN pip install \
        'h5py==3.13.*' \
        'hictkpy[all]==1.3.*' \
        'matplotlib==3.10.*' \
        'numpy==2.*' \
        'seaborn==0.13.*'

RUN python3 -c 'import h5py,hictkpy,matplotlib,numpy,seaborn'

CMD ["/bin/bash"]
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/nchg-nf'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/nchg-nf'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/nchg-nf'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-py-utils}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
