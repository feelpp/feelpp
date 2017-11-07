#!/bin/bash
# Script to push generated singularity images into the laboratory
# gitlab server. Gitlab LFS is used to store images (binaries)
# Script arguments must be a list of docker container names, for example:
# ./singularity-push.sh feelpp/feelpp-crb:latest feelpp/feelpp-toolboxes:latest

set -euo pipefail

source $(dirname $0)/common.sh

ROOT_DIR=`pwd`
IMAGES_ROOT=./singularity/images
REPO_NAME=feelpp-singularity-images.git
REPO=git@gitlab.math.unistra.fr:feelpp/feelpp-singularity-images.git
REPO_MAX_IMAGES=20 #( x10 GB )

# Remove old local images (NOT LFS).
if [ ! -d ${IMAGES_ROOT} ]; then
    echo "${IMAGES_ROOT} does not exist"
    pwd
    exit 1
fi

# Clone/update gitlab repository.
echo "-- Clone Feel++ singularity images repository"
if [ -d ${REPO_NAME} ]; then
    cd ${REPO_NAME}
    git lfs pull
    cd ${ROOT_DIR}
else
    git lfs clone --depth=1 ${REPO} ${REPO_NAME}
fi

# Copy (overide) generated images in the repository.
echo "-- Prepare Feel++ singularity images"
find ${IMAGES_ROOT} -name "*.img" \
-exec echo "- image:" {} \; \
-exec cp -r {} ${REPO_NAME} \;

cd ${REPO_NAME}
echo "-- Push Feel++ singularity images"
git add *.img
git commit -am "[buildkite] Deploy new singularity images"
git lfs push origin master
git push
cd ${ROOT_DIR}

echo "-- Remove temporary Feel++ singularity images"
find ${IMAGES_ROOT} -name "*.img" -exec  rm -r {} \;

if [ `ls -l ${REPO_NAME}/*.img | wc -l` -ge ${REPO_MAX_IMAGES} ]; then
    echo "WARNING: Maximum images limit reached ! => ${REPO}"
    echo "INFO: You should remove some images or increase this 
    REPO_MAX_IMAGES number in the script."
fi
