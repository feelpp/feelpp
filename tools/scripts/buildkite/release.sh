#!/bin/bash

set -euo pipefail

#if [ ! -z "$DOCKER_LOGIN" -a ! -z "$DOCKER_PASSWORD"]; then
if [ ! -z ${DOCKER_LOGIN+x} -a ! -z ${DOCKER_PASSWORD+x}   ]; then
    docker login --username="${DOCKER_LOGIN}" --password="${DOCKER_PASSWORD}";
fi

CONTAINERS=${*:-feelpp-libs}

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
    echo "--- Pushing Container feelpp/${container}"

    tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
        
    echo "--- Pushing feelpp/${container}:$tag"
    docker push "feelpp/${container}:$tag"

    if [ "${BUILDKITE_BRANCH}" = "develop" ]; then
        extratags=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(extratags_from_target $TARGET))
        for aliastag in ${extratags[@]} ; do
            echo "--- Pushing feelpp/${container}:$aliastag"
            docker tag "feelpp/${container}:$tag" "feelpp/${container}:$aliastag"
            docker push "feelpp/${container}:$aliastag"
        done
    fi
done

echo -e "\033[33;32m--- All tags released!\033[0m"
