#!/usr/bin/env bash
# -*- mode: sh; indent-tabs-mode: nil; sh-basic-offset: 4 -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

set -euo pipefail

version=${1:?Expected preCICE version}
checksum=${2:?Expected preCICE package SHA-256}
temporaryDir=$(mktemp -d)
trap 'rm -rf "$temporaryDir"' EXIT
packagePath="$temporaryDir/precice.deb"
curl --fail --show-error --silent --location --retry 3 \
    "https://github.com/precice/precice/releases/download/v${version}/libprecice3_${version}_noble.deb" \
    -o "$packagePath"
printf '%s  %s\n' "$checksum" "$packagePath" | sha256sum --check --strict
[[ $(dpkg-deb -f "$packagePath" Package) == libprecice3 ]]
[[ $(dpkg-deb -f "$packagePath" Version) == "$version" ]]
[[ $(dpkg-deb -f "$packagePath" Architecture) == amd64 ]]
chmod 0755 "$temporaryDir"
chmod 0644 "$packagePath"

apt-get update
apt-get install -y --no-install-recommends "$packagePath"
ldconfig
[[ $(pkg-config --modversion libprecice) == "$version" ]]
precice-version
printf 'preCICE\t%s\t%s\n' "$version" "$checksum" \
    > /usr/local/share/openfoam-precice/precice-package.tsv
