#!/usr/bin/env bash
# -*- mode: sh; indent-tabs-mode: nil; sh-basic-offset: 4 -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

set -euo pipefail

packageVersion=${1:?Expected OpenFOAM package version}
source /etc/os-release
if [[ ${ID} != ubuntu || ${VERSION_CODENAME} != noble || $(dpkg --print-architecture) != amd64 ]]; then
    echo 'This recipe currently supports Ubuntu Noble amd64 only.' >&2
    exit 1
fi

temporaryDir=$(mktemp -d)
trap 'rm -rf "$temporaryDir"' EXIT
mkdir -m 0700 "$temporaryDir/gnupg"
curl --fail --show-error --silent --location --retry 3 \
    https://dl.openfoam.com/pubkey.gpg -o "$temporaryDir/openfoam.asc"
fingerprint=$(gpg --batch --no-options --homedir "$temporaryDir/gnupg" \
    --with-colons --show-keys "$temporaryDir/openfoam.asc" | awk -F: '$1 == "fpr" { print $10; exit }')
if [[ $fingerprint != DC93C096174122E256DA24063386DD74948D208F ]]; then
    echo "Unexpected OpenFOAM signing key: $fingerprint" >&2
    exit 1
fi
gpg --batch --no-options --homedir "$temporaryDir/gnupg" --dearmor \
    --output /usr/share/keyrings/openfoam.gpg "$temporaryDir/openfoam.asc"
cat > /etc/apt/sources.list.d/openfoam.sources <<'EOF'
Types: deb
URIs: https://dl.openfoam.com/repos/deb
Suites: noble
Components: main
Architectures: amd64
Signed-By: /usr/share/keyrings/openfoam.gpg
EOF

apt-get update
apt-get install -y --no-install-recommends \
    "openfoam2606-default=$packageVersion" "openfoam2606-dev=$packageVersion"
test -f /usr/lib/openfoam/openfoam2606/etc/bashrc
mkdir -p /usr/local/share/openfoam-precice
dpkg-query -W -f='${Package}\t${Version}\t${Architecture}\n' 'openfoam2606*' \
    > /usr/local/share/openfoam-precice/openfoam-packages.tsv
cat /usr/local/share/openfoam-precice/openfoam-packages.tsv
