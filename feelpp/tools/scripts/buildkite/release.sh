#!/bin/bash

set -euo pipefail

source $(dirname $0)/common.sh

#set -x
#if [ -v DOCKER_PASSWORD -a -v DOCKER_LOGIN ]; then
#    docker login --username="${DOCKER_LOGIN}" --password="${DOCKER_PASSWORD}";
#fi
#echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin
    

build="$(basename "$0")"
VERSION=${FEELPP_VERSION}
if [ -z ${TARGET:-""} ]; then
    TARGET=ubuntu:17.04
fi
latest=false
noop=false

usage() {
    echo >&2 "usage: $build "
    echo >&2 "              [-t|--target target host os, default: $TARGET]"
    echo >&2 "              [-b|--branch project git branch, default: $BRANCH]"
    echo >&2 "              [--version project version, default: $VERSION]"
    echo >&2 "              [--latest true|false, default=$latest]"
    echo >&2 "   ie: $build -t ubuntu:16.10 -b develop --latest true feelpp-libs feelpp-base"
    exit 1
}

fromos=
fromtag=

tag=
while [ -n "$1" ]; do
    case "$1" in
        -t|--target) TARGET="$2" ; shift 2 ;;
        -b|--branch) BRANCH="$2" ; shift 2 ;;
        --version) VERSION="$2" ; shift 2 ;;
        --latest) latest="$2" ; shift 2 ;;
        --noop) noop="$2" ; shift 2 ;;
        -h|--help) usage ;;
        --) shift ; break ;;
    esac
done

BRANCHTAG=$(echo "${BRANCH}" | sed -e 's/\//-/g')
CONTAINERS=${*:-feelpp-libs}
echo $CONTAINERS

echo $CR_PAT | docker login ghcr.io -u $CR_LOGIN --password-stdin
echo "-- registered to ghcr.io..."

for container in ${CONTAINERS}; do
    echo "--- Pushing Container ghcr.io/feelpp/${container}"

    #tag=$(echo "${BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
    #tag=$(cut -d- -f 2- <<< $(tag_from_target $TARGET))
    tag=$(tag_from_target $TARGET $BRANCHTAG $VERSION)
        
    echo "--- Pushing ghcr.io/feelpp/${container}:$tag"

    if [ "$noop" = "false" ]; then
        docker push "ghcr.io/feelpp/${container}:$tag";
        docker push "ghcr.io/feelpp/${container}:$tag";
        docker push "ghcr.io/feelpp/${container}:$tag";
    else
        echo "docker push \"ghcr.io/feelpp/${container}:$tag\"";
    fi

    
    # if [ "${latest}" = "true" ]; then
    #     echo "--- Tagging ${container}:${tag}"
    #     echo "Tagging feelpp/${container}:$tag as feelpp/${container}:latest"
    #     if [ "$noop" = "false" ]; then
    #         docker tag "feelpp/${container}:$tag" "feelpp/${container}:latest";
    #         docker push  "feelpp/${container}:latest";
    #     else
    #         echo "docker tag \"feelpp/${container}:$tag\" \"feelpp/${container}:latest\"";
    #         echo "docker push  \"feelpp/${container}:latest\"";
    #     fi
    # fi
    #echo "extra"
    #if [ "${BRANCH}" = "develop" -o  "${BRANCH}" = "master" ]; then
#    if 
        #echo "extra"
        #extratags=$(echo "${BRANCH}" | sed -e 's/\//-/g')-$(cut -d- -f 2- <<< $(extratags_from_target $TARGET))
        #extratags=$(cut -d- -f 2- <<< $(extratags_from_target $TARGET))
        echo $TARGET
        echo $BRANCHTAG
        echo $VERSION
        extratags=$(echo $(extratags_from_target $TARGET $BRANCHTAG $VERSION) | tr ' ' '\n' | sort | uniq | xargs)
        echo $extratags
        for aliastag in ${extratags[@]} ; do
            echo "--- Pushing ghcr.io/feelpp/${container}:$aliastag"
            if [ "$noop" = "false" ]; then
                docker tag "ghcr.io/feelpp/${container}:$tag" "ghcr.io/feelpp/${container}:$aliastag";
                docker push "ghcr.io/feelpp/${container}:$aliastag";
                docker push "ghcr.io/feelpp/${container}:$aliastag";
                docker push "ghcr.io/feelpp/${container}:$aliastag";
            else
                echo "docker tag \"ghcr.io/feelpp/${container}:$tag\" \"ghcr.io/feelpp/${container}:$aliastag\"";
                echo "docker push \"ghcr.io/feelpp/${container}:$aliastag\"";
            fi
        done
#    fi
done

echo -e "\033[33;32m--- All tags released!\033[0m"
