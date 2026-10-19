#!/usr/bin/env bash
# -*- mode: sh; indent-tabs-mode: nil; sh-basic-offset: 4 -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

# Install the distribution's verified binary release, or build the same release on arm64.
set -euo pipefail

version=${1:?Expected preCICE version}
packageUrl=${2:?Expected preCICE package URL}
packageChecksum=${3:?Expected preCICE package SHA-256}
sourceUrl=${4:?Expected preCICE source URL}
sourceChecksum=${5:?Expected preCICE source SHA-256}
architecture=$(dpkg --print-architecture)
temporaryDir=$(mktemp -d)
trap 'rm -rf "$temporaryDir"' EXIT

downloadVerified()
{
    local url=$1 checksum=$2 destination=$3
    curl --fail --show-error --silent --location --retry 3 "$url" -o "$destination"
    printf '%s  %s\n' "$checksum" "$destination" | sha256sum --check --strict
}

case "$architecture" in
    amd64)
        packagePath="$temporaryDir/precice.deb"
        downloadVerified "$packageUrl" "$packageChecksum" "$packagePath"
        [[ $(dpkg-deb -f "$packagePath" Package) == libprecice3 ]]
        [[ $(dpkg-deb -f "$packagePath" Version) == "$version" ]]
        [[ $(dpkg-deb -f "$packagePath" Architecture) == amd64 ]]
        chmod 0755 "$temporaryDir"
        chmod 0644 "$packagePath"
        apt-get update
        apt-get install -y --no-install-recommends "$packagePath"
        installedChecksum=$packageChecksum
        ;;
    arm64)
        archivePath="$temporaryDir/precice.tar.gz"
        downloadVerified "$sourceUrl" "$sourceChecksum" "$archivePath"
        mkdir "$temporaryDir/source"
        tar -xzf "$archivePath" --strip-components=1 -C "$temporaryDir/source"
        cmake -S "$temporaryDir/source" -B "$temporaryDir/build" -G Ninja \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local \
            -DCMAKE_CXX_COMPILER=mpicxx -DMPI_CXX_COMPILER=mpicxx \
            -DPRECICE_FEATURE_MPI_COMMUNICATION=ON \
            -DPRECICE_FEATURE_PETSC_MAPPING=ON \
            -DPRECICE_FEATURE_PYTHON_ACTIONS=ON \
            -DPRECICE_BUILD_TOOLS=ON -DBUILD_TESTING=OFF
        cmake --build "$temporaryDir/build" --parallel "${PRECICE_BUILD_JOBS:-4}"
        cmake --install "$temporaryDir/build"
        installedChecksum=$sourceChecksum
        ;;
    *)
        printf 'Unsupported preCICE architecture: %s\n' "$architecture" >&2
        exit 1
        ;;
esac

ldconfig
[[ $(pkg-config --modversion libprecice) == "$version" ]]
precice-version
mkdir -p /usr/local/share/feelpp/image
printf 'preCICE\t%s\t%s\t%s\n' "$version" "$architecture" "$installedChecksum" \
    > /usr/local/share/feelpp/image/precice.tsv
apt-get clean
rm -rf /var/lib/apt/lists/*
