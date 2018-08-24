#!/bin/bash

set -eo pipefail

component=${1:-base}

source $(dirname $0)/common.sh

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

#tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
#tag=$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
BRANCHTAG=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')
tag=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
echo "--- Building feelpp-${component}:${tag}"

if [ "${component}" = "base" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-libs:${tag}" > docker/feelpp-${component}/dockerfile.tmp
elif [ "${component}" = "pyfeelpp" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-libs:${tag}" > docker/feelpp-${component}/dockerfile.tmp
elif [ "${component}" = "toolboxes" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-base:${tag}" > docker/feelpp-${component}/dockerfile.tmp
elif [ "${component}" = "pyfeelpp-toolboxes" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-toolboxes:${tag}" > docker/feelpp-${component}/dockerfile.tmp    
elif [ "${component}" = "mor" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-toolboxes:${tag}" > docker/feelpp-${component}/dockerfile.tmp
elif [ "${component}" = "pyfeelpp-mor" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-mor:${tag}" > docker/feelpp-${component}/dockerfile.tmp
else
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-toolboxes:${tag}" > docker/feelpp-${component}/dockerfile.tmp
fi    
docker build \
       --tag=feelpp/feelpp-${component}:${tag} \
       --build-arg=BUILD_JOBS=${JOBS}\
       --build-arg=BRANCH=${BUILDKITE_BRANCH}\
       --build-arg=CXX="${CXX}"\
       --build-arg=CC="${CC}" \
       --no-cache=true \
       -f docker/feelpp-${component}/dockerfile.tmp \
       docker/feelpp-${component}


echo "--- Tagging feelpp-${component}:${tag}"
extratags=$(extratags_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
# add extra tags
for tagalias in ${extratags[@]}; do
    echo "Tagging feelpp/feelpp-${component}:$tag as feelpp/feelpp-${component}:$tagalias"
    docker tag "feelpp/feelpp-${component}:$tag" "feelpp/feelpp-${component}:$tagalias"
done

    


