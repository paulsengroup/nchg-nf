# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

ARG CONTAINER_VERSION

FROM paulsengroup/hictk:${CONTAINER_VERSION} AS base

ARG CONTAINER_TITLE
ARG CONTAINER_VERSION

ARG PIP_NO_CACHE_DIR=0

RUN apt-get update \
&& apt-get install -y procps \
&& rm -rf /var/lib/apt/lists/*


RUN hictk --help

CMD ["/usr/local/bin/hictk"]
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/nchg-nf'
LABEL org.opencontainers.image.documentation='https://github.com/paulsengroup/nchg-nf'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/nchg-nf'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hictk}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
