#!/bin/sh

containers="$1"

set -euo pipefail

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

cd ./docker/singularity

for cont in ${containers}
do
    echo "--- generate singularity image for ${cont} docker container"
    ./generate_bootstrap.sh "${cont}"
    ./generate_image.sh "${cont}"
done
