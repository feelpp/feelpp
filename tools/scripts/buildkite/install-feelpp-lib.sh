#!/bin/bash

set -eo pipefail

tag_from_target() {
    splitfrom=(`echo "$TARGET" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}

    tools/scripts/buildkite/list.sh | grep "${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        image=${tokens[0]}
        printf "%s" "$image" 
    done
}
extratags_from_target() {
    splitfrom=(`echo "$TARGET" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}

    tools/scripts/buildkite/list.sh | grep "${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        extratags=${tokens[@]:5}
        printf "%s" "${extratags}" 
    done
}

echo '--- clone/pull feelpp/docker'
if [ -d docker ]; then (cd docker; git pull) else git clone --depth=1 https://github.com/feelpp/docker; fi

echo '--- building feelpp-libs'
(cd docker/feelpp-libs && bash mkimg.sh -f ${TARGET} --jobs ${JOBS} --branch ${BUILDKITE_BRANCH} --cxx "${CXX}" --cc "${CC}")

extratags=`extratags_from_target $TARGET`
# add extra tags
for tagalias in ${extratags[@]}; do
    echo "Tagging feelpp/feelpp-libs:$tag as feelpp/feelpp-libs:$tagalias"
    docker tag "feelpp/feelpp-libs:$tag" "feelpp/feelpp-libs:$tagalias"
done

