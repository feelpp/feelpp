#!/bin/bash

set -euo pipefail

#if [ ! -z "$DOCKER_LOGIN" -a ! -z "$DOCKER_PASSWORD"]; then
#    docker login --username="${DOCKER_LOGIN}" --password="${DOCKER_PASSWORD}";
#fi

CONTAINERS=${*:-libs}

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

for container in ${CONTAINERS}; do
    echo "--- Pushing Container feelpp/feelpp-${container}"

    tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
        
    echo "--- Pushing feelpp/feelpp-${container}:$tag"
    docker push "feelpp/feelpp-${container}:$tag"

    for aliastag in ${extratags[@]} ; do
        echo "--- Pushing feelpp/feelpp-${container}:$aliastag"
        docker tag "feelpp/feelpp-${container}:$tag" "feelpp/feelpp-${container}:$aliastag"
        docker push "feelpp/feelpp-${container}:$aliastag"
    done
done

echo -e "\033[33;32m--- All tags released!\033[0m"
