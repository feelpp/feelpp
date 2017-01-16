#!/bin/bash

set -eo pipefail

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

echo '--- building feelpp-libs'
cd docker/feelpp-libs && bash mkimg.sh -f ${TARGET} --jobs ${JOBS} --branch ${BUILDKITE_BRANCH} --cxx "${CXX}" --cc "${CC}"

