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

echo '--- building feelpp-base'
tag=`tag_from_target $TARGET`

cd docker/feelpp-base \
    && docker build \
              --tag=feelpp/feelpp-base:${tag} \
              --build-arg=BUILD_JOBS=${JOBS}\
              --build-arg=BRANCH=${BUILDKITE_BRANCH}\
              --build-arg=CXX="${CXX}"\
              --build-arg=CC="${CC}" \
              . 
extratags=`extratags_from_target $TARGET`
# add extra tags
for tagalias in ${extratags[@]}; do
    echo "Tagging feelpp/feelpp-base:$tag as feelpp/feelpp-base:$tagalias"
    docker tag "feelpp/feelpp-base:$tag" "feelpp/feelpp-base:$tagalias"
done

    


