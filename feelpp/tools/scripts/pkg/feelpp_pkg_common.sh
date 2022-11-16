#! /bin/bash

set -eo pipefail

BUILDKITE_AGENT_NAME=${BUILDKITE_AGENT_NAME:-default}
BUILDKITE_BRANCH=${BUILDKITE_BRANCH:-develop}
BRANCH=${BRANCH:-${BUILDKITE_BRANCH}}

# default values
CHANNEL=${CHANNEL:-latest}
if [ "$BUILDKITE_BRANCH" = "develop" -o  "$BRANCH" = "develop" ]; then
    CHANNEL=latest
fi
if [ "$BUILDKITE_BRANCH" = "master" -o  "$BRANCH" = "master" ]; then
    CHANNEL=stable
fi 
DIST=${DIST:-focal}
if [ "$DIST" = "bionic" -o "$DIST" = "eoan" -o "$DIST" = "focal" -o "$DIST" = "jammy" -o "$DIST" = "kinetic" ]; then
   FLAVOR=ubuntu
elif [ "$DIST" = "buster" -o "$DIST" = "bullseye" -o "$DIST" = "bookworm" -o "$DIST" = "testing" -o "$DIST" = "sid" ]; then
    FLAVOR=debian
fi



COMPONENT=${COMPONENT:-feelpp}