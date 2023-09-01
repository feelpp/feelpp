#!/bin/bash

set -eo pipefail
set -x
component=${1:-base}

source $(dirname $0)/common.sh

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

#tag=$(echo "${GITHUB_REF_NAME}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
#tag=$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
BRANCHTAG=$(echo "${GITHUB_REF_NAME}" | sed -e 's/\//-/g')
tag_compiler=$(echo ${CC} | sed -e 's/-//g')
if test "$tag_compiler" != "${tag_compiler%gcc*}"; then
    tag=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)-${tag_compiler}
else
    tag=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
fi
tagos=$(tag_from_os $TARGET $BRANCHTAG $FEELPP_VERSION)
if test "$BUILDKITE_PIPELINE_SLUG" = "feelpp-debug"; then
    tag=${tag}-debug
fi
image="feelpp-${component}"
if [ "${component}" = "feelpp" ] ; then
#    tag=$(tag_from_os $TARGET $BRANCHTAG $FEELPP_VERSION)
    image="feelpp"
fi
if [ "${component}" = "feelpp-python" ] ; then
#    tag=$(tag_from_os $TARGET $BRANCHTAG $FEELPP_VERSION)
    image="feelpp-python"
fi
echo "--- Building ${image}:${tag}"


if [ "${component}" = "feelpp" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-env:${tagos}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "toolboxes" -o "${component}" = "testsuite" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp:${tag}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "mor" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-toolboxes:${tag}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "feelpp-python" -o "${component}" = "python" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-mor:${tag}" > docker/${image}/dockerfile.tmp
else
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-toolboxes:${tag}" > docker/${image}/dockerfile.tmp
fi

if [ "${component}" = "feelpp" ] ; then
    CTEST_FLAGS="-R feelpp_qs_ -T test --no-compress-output"
elif [ "${component}" = "toolboxes" ] ; then
    CTEST_FLAGS="-R feelpp_toolbox_ -T test --no-compress-output --output-on-failure"
elif [ "${component}" = "testsuite" ] ; then
    CTEST_FLAGS="-R feelpp_test_ -T test --no-compress-output --output-on-failure"
elif [ "${component}" = "feelpp-python" -o "${component}" = "python" ] ; then
    CTEST_FLAGS="-R feelpp -T test --no-compress-output --output-on-failure"
else
    CTEST_FLAGS="-T test --no-compress-output --output-on-failure"
fi

docker build \
       --pull \
       --tag=ghcr.io/feelpp/${image}:${tag} \
       --build-arg=BUILD_JOBS=${JOBS}\
       --build-arg=BRANCH=${GITHUB_REF_NAME}\
       --build-arg=CXX="${CXX}" \
       --build-arg=CC="${CC}" \
       --build-arg=CONFIGURE_FLAGS="${CONFIGURE_FLAGS}" \
       --build-arg=CMAKE_FLAGS="${CMAKE_FLAGS}" \
       --build-arg=CTEST_FLAGS="${CTEST_FLAGS}" \
       --no-cache=true \
       -f docker/${image}/dockerfile.tmp \
       docker/${image}


echo "--- Tagging ghcr.io/feelpp/${image}:${tag}"
extratags=$(extratags_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
# add extra tags
for tagalias in ${extratags[@]}; do
    echo "Tagging ghcr.io/feelpp/${image}:$tag as ghcr.io/feelpp/${image}:$tagalias"
    docker tag "ghcr.io/feelpp/${image}:$tag" "ghcr.io/feelpp/${image}:$tagalias"
done
source $(dirname $0)/release.sh  -- ${image}
