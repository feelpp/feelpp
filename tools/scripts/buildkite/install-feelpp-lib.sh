#!/bin/bash

set -eo pipefail

tag_from_target() {
    splitfrom=(`echo "$TARGET" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}

    tools/scripts/buildkite/list.sh | grep "${BUILDKITE_BRANCH}-${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        image=${tokens[0]}
        printf "%s" "${image}"
    done
}
extratags_from_target() {
    splitfrom=(`echo "$TARGET" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}

    tools/scripts/buildkite/list.sh | grep "${BUILDKITE_BRANCH}-${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        extratags=${tokens[@]:5}
        printf "%s" "${extratags}" 
    done
}

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

#tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
#tag=$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
tag=$(tag_from_target $TARGET)

echo "--- building feelpp-libs:${tag}"
(cd docker/feelpp-libs && bash mkimg.sh -f ${TARGET} --jobs ${JOBS} --branch ${BUILDKITE_BRANCH} --cxx "${CXX}" --cc "${CC}")

if [ "${BUILDKITE_BRANCH}" = "develop" ]; then
    echo "--- Tagging feelpp-libs:${tag}"
    extratags=`extratags_from_target $TARGET`
    # add extra tags
    for tagalias in ${extratags[@]}; do
        echo "Tagging feelpp/feelpp-libs:$tag as feelpp/feelpp-libs:$tagalias"
        docker tag "feelpp/feelpp-libs:$tag" "feelpp/feelpp-libs:$tagalias"
    done
fi
