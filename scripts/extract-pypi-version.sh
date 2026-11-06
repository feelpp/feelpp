#!/bin/bash
# Extract Feel++ version from feelpp.version.cmake and convert to PEP 440
# Usage: ./scripts/extract-pypi-version.sh [path/to/feelpp.version.cmake]
#
# PEP 440 mapping:
#   0.111.0-preview.12  ->  0.111.0.dev12
#   0.111.0-rc.1        ->  0.111.0rc1
#   0.111.0-alpha.3     ->  0.111.0a3
#   0.111.0-beta.2      ->  0.111.0b2
#   0.111.0             ->  0.111.0

set -euo pipefail

VERSION_FILE="${1:-feelpp.version.cmake}"

if [ ! -f "$VERSION_FILE" ]; then
    echo "Error: $VERSION_FILE not found" >&2
    exit 1
fi

MAJOR=$(grep 'set(FEELPP_VERSION_MAJOR' "$VERSION_FILE" | sed 's/.*"\(.*\)".*/\1/')
MINOR=$(grep 'set(FEELPP_VERSION_MINOR' "$VERSION_FILE" | sed 's/.*"\(.*\)".*/\1/')
MICRO=$(grep 'set(FEELPP_VERSION_MICRO' "$VERSION_FILE" | sed 's/.*"\(.*\)".*/\1/')
PRE=$(grep 'set(FEELPP_VERSION_PRERELEASE' "$VERSION_FILE" | sed 's/.*"\(.*\)".*/\1/')

VERSION="${MAJOR}.${MINOR}.${MICRO}"

case "$PRE" in
    *preview*)
        NUM=$(echo "$PRE" | grep -oP '\d+$')
        VERSION="${VERSION}.dev${NUM}"
        ;;
    *rc*)
        NUM=$(echo "$PRE" | grep -oP '\d+$')
        VERSION="${VERSION}rc${NUM}"
        ;;
    *alpha*)
        NUM=$(echo "$PRE" | grep -oP '\d+$')
        VERSION="${VERSION}a${NUM}"
        ;;
    *beta*)
        NUM=$(echo "$PRE" | grep -oP '\d+$')
        VERSION="${VERSION}b${NUM}"
        ;;
esac

echo "$VERSION"
