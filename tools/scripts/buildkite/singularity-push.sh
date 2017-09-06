#!/bin/bash

set -euo pipefail

source $(dirname $0)/common.sh

REPO_NAME=feelpp-singularity-images.git

# Remove old images if exist.
rm -r ./docker.git/singularity/images/*

# Clone/update gitlab repository.
if [ -d ${REPO_NAME} ]; then
	cd ${REPO_NAME}
	git lfs pull
else
	git lfs clone --depth=1 git@gitlab.math.unistra.fr:dolle/feelpp-singularity-images.git ${REPO_NAME}
fi

# Copy (overide) generated images in the repository.
find ./docker.git/singularity/images/ -name "*.img" \
-exec echo "--- push singularity image:" {} \; \
-exec cp -r {} ${REPO_NAME} \;

cd ${REPO_NAME}
git add *.img
git commit -am "[buildkite] Deploy new singularity images"
git push
