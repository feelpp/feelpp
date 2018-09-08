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

image="feelpp-${component}"
if [ "${component}" = "feelpp" ] ; then
    image="feelpp"
fi

if [ "${component}" = "feelpp" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "feelpp/feelpp-env:${tag}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "toolboxes" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "feelpp/feelpp:${tag}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "mor" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "feelpp/feelpp-toolboxes:${tag}" > docker/${image}/dockerfile.tmp
else
    dockerfile_from "docker/${image}/Dockerfile.template" "feelpp/feelpp-toolboxes:${tag}" > docker/${image}/dockerfile.tmp
fi    
docker build \
       --tag=feelpp/${image}:${tag} \
       --build-arg=BUILD_JOBS=${JOBS}\
       --build-arg=BRANCH=${BUILDKITE_BRANCH}\
       --build-arg=CXX="${CXX}"\
       --build-arg=CC="${CC}" \
       --no-cache=true \
       -f docker/${image}/dockerfile.tmp \
       docker/${image}


echo "--- Tagging ${image}:${tag}"
extratags=$(extratags_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
# add extra tags
for tagalias in ${extratags[@]}; do
    echo "Tagging feelpp/${image}:$tag as feelpp/${image}:$tagalias"
    docker tag "feelpp/${image}:$tag" "feelpp/${image}:$tagalias"
done

    


