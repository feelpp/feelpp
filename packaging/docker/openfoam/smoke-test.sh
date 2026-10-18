#!/usr/bin/env bash
# -*- mode: sh; indent-tabs-mode: nil; sh-basic-offset: 4 -*-
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-License-Identifier: LGPL-3.0-or-later

set -euo pipefail
if (( $# > 1 )) || { (( $# == 1 )) && [[ $1 != --mpi ]]; }; then
    echo 'Usage: openfoam-smoke-test [--mpi]' >&2
    exit 2
fi

foamApi=$(foamEtcFile -show-api)
[[ $foamApi == 2606 ]]
printf 'OpenFOAM API: %s\n' "$foamApi"
precice-version
ompi_info --version
cat /usr/local/share/openfoam-precice/*.tsv

temporaryDir=$(mktemp -d)
trap 'rm -rf "$temporaryDir"' EXIT
adapterPath="${PRECICE_OPENFOAM_TARGET_DIR}/libpreciceAdapterFunctionObject.so"
for library in "$adapterPath" "$(command -v laplacianFoam)" "$(command -v openfoam-mpi-probe)"; do
    linkStatus=0
    ldd -r "$library" > "$temporaryDir/ldd.log" 2>&1 || linkStatus=$?
    cat "$temporaryDir/ldd.log"
    if (( linkStatus != 0 )) || grep -E 'not found|undefined symbol' "$temporaryDir/ldd.log"; then
        exit 1
    fi
done
foamHasLibrary -verbose precice libpreciceAdapterFunctionObject
laplacianFoam -help
pimpleFoam -help

if [[ ${1:-} == --mpi ]]; then
    mpiexec -n 2 --bind-to core --report-bindings openfoam-mpi-probe
fi
