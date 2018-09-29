#!/bin/bash

set -eo pipefail

component=${1:-base}

source $(dirname $0)/common.sh

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

#tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
#tag=$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
BRANCHTAG=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')
tag_compiler=$(echo ${CC} | sed -e 's/-//g')
if test "$tag_compiler" != "${tag_compiler%gcc*}"; then
    tag=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)-${tag_compiler}
else
    tag=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
fi
image="feelpp-${component}"
if [ "${component}" = "feelpp" ] ; then
#    tag=$(tag_from_os $TARGET $BRANCHTAG $FEELPP_VERSION)
    image="feelpp"
fi
echo "--- Building ${image}:${tag}"


if [ "${component}" = "feelpp" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "feelpp/feelpp-env:latest" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "toolboxes" -o "${component}" = "testsuite" ] ; then
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

    


