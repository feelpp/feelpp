#!/bin/bash

set -euo pipefail

source $(dirname $0)/common.sh

IMAGES_ROOT=./docker/singularity/images
REPO_NAME=feelpp-singularity-images.git

# Remove old local images (NOT LFS).
if [ -d ${IMAGES_ROOT} ]; then
    find ${IMAGES_ROOT} -name "*.img" -exec  rm -r {} \;
else
    echo "./docker.git/singularity/images does not exist"
    pwd
    exit 1
fi

# Clone/update gitlab repository.
if [ -d ${REPO_NAME} ]; then
	cd ${REPO_NAME}
	git lfs pull
else
	git lfs clone --depth=1 git@gitlab.math.unistra.fr:dolle/feelpp-singularity-images.git ${REPO_NAME}
fi

# Copy (overide) generated images in the repository.
find ${IMAGES_ROOT} -name "*.img" \
-exec echo "--- push singularity image:" {} \; \
-exec cp -r {} ${REPO_NAME} \;

cd ${REPO_NAME}
git add *.img
git commit -am "[buildkite] Deploy new singularity images"
git push
