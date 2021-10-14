#! /bin/bash

set -eo pipefail

BUILDKITE_AGENT_NAME=${BUILDKITE_AGENT_NAME:-default}
# default values
CHANNEL=latest
if [ "$BUILDKITE_BRANCH" = "develop" -o  "$BRANCH" = "develop" ]; then
    CHANNEL=latest
fi
if [ "$BUILDKITE_BRANCH" = "master" -o  "$BRANCH" = "master" ]; then
    CHANNEL=stable
fi 
DIST=${DIST:-focal}
if [ "$DIST" = "bionic" -o "$DIST" = "eoan" -o "$DIST" = "focal"  -o "$DIST" = "groovy"  -o "$DIST" = "hirsute"  ]; then
   FLAVOR=ubuntu
elif [ "$DIST" = "buster" -o "$DIST" = "bullseye" -o "$DIST" = "sid" ]; then
    FLAVOR=debian
fi



COMPONENT=${COMPONENT:-feelpp}