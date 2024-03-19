#!/bin/bash
# detect_platform.sh

ARCH=$(uname -m)
if [ "$ARCH" = "x86_64" ]; then
    export PLATFORM_TYPE="x86"
elif [ "$ARCH" = "aarch64" ]; then
    export PLATFORM_TYPE="arm64"
else
    export PLATFORM_TYPE="unknown"
fi

echo "Detected platform: $PLATFORM_TYPE"
# Optionally write the platform type to a file or perform other actions
