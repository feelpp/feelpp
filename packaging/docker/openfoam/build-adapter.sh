#!/usr/bin/env bash
# -*- mode: sh; indent-tabs-mode: nil; sh-basic-offset: 4 -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

set -euo pipefail

version=${1:?Expected adapter version}
checksum=${2:?Expected adapter source SHA-256}
buildJobs=${3:?Expected parallel build count}
temporaryDir=$(mktemp -d)
trap 'rm -rf "$temporaryDir"' EXIT
curl --fail --show-error --silent --location --retry 3 \
    "https://github.com/precice/openfoam-adapter/archive/refs/tags/v${version}.tar.gz" \
    -o "$temporaryDir/adapter.tar.gz"
printf '%s  %s\n' "$checksum" "$temporaryDir/adapter.tar.gz" | sha256sum --check --strict
tar -xzf "$temporaryDir/adapter.tar.gz" -C "$temporaryDir"

# OpenFOAM's setup probes optional libraries with nonzero return statuses and
# expects unset variables. Restore strict checks after loading its environment.
set +eu
source "${FOAM_BASHRC}" '' || exit $?
set -eu
[[ ${WM_PROJECT_VERSION} == v2606 || ${WM_PROJECT_VERSION} == 2606 ]]
export WM_NCOMPPROCS="$buildJobs"
mkdir -p "${PRECICE_OPENFOAM_TARGET_DIR}"
cd "$temporaryDir/openfoam-adapter-$version"
logDir=/usr/local/share/openfoam-precice/build-logs
mkdir -p "$logDir"
buildStatus=0
./Allwmake || buildStatus=$?
for logFile in Allwmake.log wmake.log ldd.log; do
    if [[ -f $logFile ]]; then
        cp "$logFile" "$logDir/"
    fi
done
if (( buildStatus != 0 )); then
    exit "$buildStatus"
fi

# Allwmake can report link errors without returning a failure status.
adapterPath="${PRECICE_OPENFOAM_TARGET_DIR}/libpreciceAdapterFunctionObject.so"
test -s "$adapterPath"
linkStatus=0
ldd -r "$adapterPath" > "$logDir/ldd-verified.log" 2>&1 || linkStatus=$?
cat "$logDir/ldd-verified.log"
if (( linkStatus != 0 )) || grep -E 'not found|undefined symbol' "$logDir/ldd-verified.log"; then
    exit 1
fi
install -m 0644 LICENSE /usr/local/share/openfoam-precice/adapter-LICENSE
printf 'openfoam-adapter\t%s\t%s\n' "$version" "$checksum" \
    > /usr/local/share/openfoam-precice/adapter-source.tsv
