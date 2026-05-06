#!/bin/bash

set -euo pipefail

# Script can be run from the Feel++ top-level directory or any other cwd.
scriptdir="$(cd "$(dirname "$0")" && pwd)"
source "$scriptdir/feelpp_pkg_common.sh"

# Ensure DIST, FLAVOR, and CHANNEL are set
if [[ -z "${DIST:-}" || -z "${FLAVOR:-}" || -z "${CHANNEL:-}" ]]; then
    echo "Error: DIST, FLAVOR, and CHANNEL must be set."
    exit 1
fi

echo "--- preparing pbuilder base for DIST=$DIST, FLAVOR=$FLAVOR, CHANNEL=$CHANNEL"
prepare_feelpp_pbuilder_base "$DIST"
