#!/bin/bash

set -eo pipefail

source $(dirname $0)/common.sh

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

#tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
#tag=$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
BRANCHTAG=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')
tag=$(tag_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
echo "--- building feelpp-libs:${tag}"

(cd docker/feelpp-libs && bash mkimg.sh -t feelpp/feelpp-libs:$tag -f ${TARGET} --jobs ${JOBS} --branch ${BUILDKITE_BRANCH} --cxx "${CXX}" --cc "${CC}")

if [ "${BUILDKITE_BRANCH}" = "develop" ]; then
    echo "--- Tagging feelpp-libs:${tag}"
    extratags=$(extratags_from_target $TARGET $BRANCHTAG $FEELPP_VERSION)
    # add extra tags
    for tagalias in ${extratags[@]}; do
        echo "Tagging feelpp/feelpp-libs:$tag as feelpp/feelpp-libs:$tagalias"
        docker tag "feelpp/feelpp-libs:$tag" "feelpp/feelpp-libs:$tagalias"
    done
fi
