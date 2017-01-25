#!/bin/bash

set -eo pipefail

component=${1:-base}

COMPONENTS=(base toolboxes crb) 

# Combines a dockerfile template with a generated FROM line
dockerfile_from() {
    local dockerfile="$1"
    local from="$2"
    printf 'FROM %s\n%s' "$from" "$(<$dockerfile)"
}

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

tag=$(echo "${BUILDKITE_BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
echo "--- Building feelpp-${component}:${tag}"

if [ "${component}" = "base" ] ; then
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-libs:${tag}" > docker/feelpp-${component}/dockerfile.tmp
else
    dockerfile_from "docker/feelpp-${component}/Dockerfile.template" "feelpp/feelpp-base:${tag}" > docker/feelpp-${component}/dockerfile.tmp
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


if [ "${BUILDKITE_BRANCH}" = "develop" ]; then
    echo "--- Tagging feelpp-${component}:${tag}"
    extratags=`extratags_from_target $TARGET`
    # add extra tags
    for tagalias in ${extratags[@]}; do
        echo "Tagging feelpp/feelpp-${component}:$tag as feelpp/feelpp-${component}:$tagalias"
        docker tag "feelpp/feelpp-${component}:$tag" "feelpp/feelpp-${component}:$tagalias"
    done
fi
    


