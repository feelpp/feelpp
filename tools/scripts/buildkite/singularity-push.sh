#!/bin/bash
# Script to push generated singularity images into the laboratory
# gitlab server. Gitlab LFS is used to store images (binaries)
# Script arguments must be a list of docker container names, for example:
# ./singularity-push.sh feelpp/feelpp-crb:latest feelpp/feelpp-toolboxes:latest

containers="$@"

set -euo pipefail

source $(dirname $0)/common.sh

BRANCHTAG=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')
FEELPP_DOCKER_TAG=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)

echo '--- clone/pull feelpp/singularity'
if [ -d singularity ]; then
    cd singularity; git pull
else
    git clone --depth=1 https://github.com/feelpp/singularity;
    cd singularity
fi

for cont in ${containers}
do
    echo "--- generate singularity image for ${cont} docker container using tag: ${FEELPP_DOCKER_TAG}"
    ./push_image_ftp_server.sh "${cont}:${FEELPP_DOCKER_TAG}"
    #./push_image_lfs_server.sh "${cont}:${FEELPP_DOCKER_TAG}"
done
