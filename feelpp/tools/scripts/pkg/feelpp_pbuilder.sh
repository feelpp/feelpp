#!/bin/bash

set -euo pipefail

# Script must be run from the Feel++ top-level directory
scriptdir="$(cd "$(dirname "$0")" && pwd)"
source "$scriptdir/feelpp_pkg_common.sh"

# Ensure DIST, FLAVOR, and CHANNEL are set
if [[ -z "${DIST:-}" || -z "${FLAVOR:-}" || -z "${CHANNEL:-}" ]]; then
    echo "Error: DIST, FLAVOR, and CHANNEL must be set."
    exit 1
fi

# Read build dependencies
if [[ ! -f "$DIST" ]]; then
    echo "Error: File '$DIST' with build dependencies not found."
    exit 1
fi
builddeps=$(cat "$DIST" | tr "\n" " ")

# Start pbuilder update
feelpp-pbuilder-dist "$DIST" login --save-after-login << EOF
echo "--- Running in pbuilder for DIST=$DIST, FLAVOR=$FLAVOR, CHANNEL=$CHANNEL"

echo "--- Updating package lists"
apt-get update -yq

echo "--- Installing base tools"
apt-get install -yq apt-transport-https ca-certificates gnupg software-properties-common wget


echo "--- Removing old Feel++ repository configurations"
rm -f /etc/apt/sources.list.d/feelpp.list

echo "--- Adding Feel++ repository"
if [ $DIST = "jammy" -o  $DIST = "focal" -o $DIST = "bookworm" ]; then
    wget -O - http://apt.feelpp.org/apt.gpg | apt-key add -
    echo 'deb [trusted=yes] http://apt.feelpp.org/ubuntu/jammy jammy latest' > /etc/apt/sources.list.d/feelpp.list 
else
    wget -O - http://apt.feelpp.org/apt.gpg 2>/dev/null | gpg --dearmor - |  tee /usr/share/keyrings/feelpp-archive-keyring.gpg >/dev/null
    echo 'deb [signed-by=/usr/share/keyrings/feelpp-archive-keyring.gpg] http://apt.feelpp.org/$FLAVOR/$DIST $DIST $CHANNEL' | tee /etc/apt/sources.list.d/feelpp.list >/dev/null
fi

echo "--- Importing Feel++ repository GPG key"
if ! grep -q "feelpp" /etc/apt/sources.list.d/feelpp.list; then
    echo "Error: Failed to add Feel++ repository."
    exit 1
fi


if [ "$DIST" = "focal" ]; then
    echo "--- Adding Kitware repository for Focal"
    wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc | gpg --dearmor -o /usr/share/keyrings/kitware-archive-keyring.gpg
    echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/$FLAVOR/ $DIST main" > /etc/apt/sources.list.d/kitware.list
fi

if [ "$DIST" = "bullseye" ]; then
    echo "--- Adding Debian Backports for Bullseye"
    echo "deb http://deb.debian.org/debian $DIST-backports main" > /etc/apt/sources.list.d/backports.list
fi

echo "--- Updating package lists after adding repositories"
apt-get update -yq

echo "--- Installing build dependencies"
echo "Build dependencies: $builddeps"
apt-get install -yq $builddeps

if [ "$DIST" = "bullseye" ]; then
    echo "--- Installing CMake from Bullseye Backports"
    apt-get install -yq -t bullseye-backports cmake
fi

echo "--- Cleaning up"
apt-get clean

EOF