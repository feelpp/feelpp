#!/bin/bash

set -euo pipefail

docker login --username="${DOCKER_LOGIN}" --password="${DOCKER_PASSWORD}"

CONTAINERS=${*:-libs}

for container in ${CONTAINERS}; do
    echo "--- Pushing Container feelpp/feelpp-${container}"
    tools/scripts/buildkite/list.sh | while read line ; do
        tokens=($line)
        image=${tokens[0]}
        extratags=${tokens[@]:5}

        echo "--- Pushing feelpp/feelpp-${container}:$image"
        docker push "feelpp/feelpp-${container}:$image"

        for aliastag in ${extratags[@]} ; do
            echo "--- Pushing feelpp/feelpp-${container}:$aliastag"
            docker tag "feelpp/feelpp-${container}:$image" "feelpp/feelpp-${container}:$aliastag"
            docker push "feelpp/feelpp-${container}:$aliastag"
        done
    done
done

echo -e "\033[33;32m--- All images released!\033[0m"
