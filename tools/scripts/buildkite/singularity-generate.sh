#!/bin/bash

containers="$@"

set -euo pipefail

source $(dirname $0)/common.sh

BRANCHTAG=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')
FEELPP_DOCKER_TAG=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

cd ./docker/singularity

for cont in ${containers}
do
    echo "--- generate singularity image for ${cont} docker container using tag: ${FEELPP_DOCKER_TAG}"
    ./generate_bootstrap.sh "${cont}"
    ./generate_image.sh "${cont}"
done
