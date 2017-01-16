#!/bin/bash

set -euo pipefail

tools/scripts/buildkite/list.sh | while read line ; do
  tokens=($line)
  image=${tokens[0]}
  extratags=${tokens[@]:5}

  echo "--- Pushing feelpp/feelpp-libs:$image"
  docker push "feelpp/feelpp-libs:$image"

  for aliastag in ${extratags[@]} ; do
    echo "--- Pushing feelpp/feelpp-libs:$aliastag"
    docker push "feelpp/feelpp-libes:$aliastag"
  done
done

echo -e "\033[33;32m--- All images released!\033[0m"
