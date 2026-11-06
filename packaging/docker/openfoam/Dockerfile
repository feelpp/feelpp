# syntax=docker/dockerfile:1
# -*- mode: Dockerfile -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

ARG BASE_IMAGE=ubuntu:24.04
FROM ${BASE_IMAGE} AS openfoam

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG DEBIAN_FRONTEND=noninteractive
ARG OPENFOAM_PACKAGE_VERSION=2606.0~rc2-1

ENV FOAM_BASHRC=/usr/lib/openfoam/openfoam2606/etc/bashrc \
    PRECICE_OPENFOAM_TARGET_DIR=/usr/local/lib/openfoam-precice

COPY install-openfoam.sh entrypoint.sh /opt/feelpp-openfoam-recipe/
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates curl gnupg build-essential pkg-config \
    && bash /opt/feelpp-openfoam-recipe/install-openfoam.sh "${OPENFOAM_PACKAGE_VERSION}" \
    && chmod 0755 /opt/feelpp-openfoam-recipe/entrypoint.sh \
    && useradd --uid 1001 --create-home --shell /bin/bash openfoam \
    && mkdir -p /work \
    && chown openfoam:openfoam /work \
    && rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.source="https://github.com/feelpp/feelpp" \
      org.opencontainers.image.vendor="University of Strasbourg" \
      org.opencontainers.image.title="OpenFOAM 2606 development environment" \
      org.opencontainers.image.version="${OPENFOAM_PACKAGE_VERSION}"

WORKDIR /work
USER openfoam
ENTRYPOINT ["/opt/feelpp-openfoam-recipe/entrypoint.sh"]
CMD ["bash"]

FROM openfoam AS openfoam-precice
USER root

ARG PRECICE_VERSION=3.4.1
ARG PRECICE_DEB_SHA256=3a36a40254acae06409d067381df1b68c309cc4ffb1b1d15590e1fa5359be888
ARG ADAPTER_VERSION=1.4.0
ARG ADAPTER_SHA256=399c89ab0fc5dfcff87f77433dddb60d4de0b2071de9b78b0e756e4079fd6a21
ARG BUILD_JOBS=4

COPY install-precice.sh build-adapter.sh smoke-test.sh mpi-probe.c /opt/feelpp-openfoam-recipe/
RUN bash /opt/feelpp-openfoam-recipe/install-precice.sh "${PRECICE_VERSION}" "${PRECICE_DEB_SHA256}" \
    && bash /opt/feelpp-openfoam-recipe/build-adapter.sh "${ADAPTER_VERSION}" "${ADAPTER_SHA256}" "${BUILD_JOBS}" \
    && mpicc -Wall -Wextra -Werror /opt/feelpp-openfoam-recipe/mpi-probe.c -o /usr/local/bin/openfoam-mpi-probe \
    && install -m 0755 /opt/feelpp-openfoam-recipe/smoke-test.sh /usr/local/bin/openfoam-smoke-test \
    && /opt/feelpp-openfoam-recipe/entrypoint.sh openfoam-smoke-test \
    && rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.title="OpenFOAM 2606 with preCICE" \
      org.opencontainers.image.description="Independent OpenFOAM participant for coupling with Feel++ and KUB" \
      org.feelpp.precice.version="${PRECICE_VERSION}" \
      org.feelpp.openfoam-adapter.version="${ADAPTER_VERSION}"

USER openfoam
