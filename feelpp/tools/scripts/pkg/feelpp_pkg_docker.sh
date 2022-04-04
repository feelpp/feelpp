#!/bin/bash

set -eo pipefail
#set -x
source $(dirname $0)/feelpp_pkg_common.sh

TOKEN="${CR_PAT}"
REPO="feelpp/docker" # format: username/repository
EVENT_TYPE="pkg-feelpp-published"

#cmake --preset feelpp
#cmake --build --preset feelpp -t dist
COMPONENT=feelpp
main_version=$(echo build/$COMPONENT/feelpp-*.tar.gz | sed  "s/build\/feelpp\/feelpp-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\1/g")
extra_version=$(echo build/$COMPONENT/feelpp-*.tar.gz | sed "s/build\/feelpp\/feelpp-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\2/g")
if [ -z $extra_version ]; then
    version=$(echo build/$COMPONENT/feelpp-*.tar.gz | sed  "s/build\/feelpp\/feelpp-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\1/g" )
else
    version=$(echo build/$COMPONENT/feelpp-*.tar.gz | sed  "s/build\/feelpp\/feelpp-\([0-9.]*\)-*\([a-z.0-9]*\).tar.gz/\1-\2/g" )
fi

echo $DIST
echo $version
echo $EVENT_TYPE
echo "{\"event_type\": \"$EVENT_TYPE\", \"client_payload\":{ \"dist\": \"$DIST\", \"version\": \"$version\", \"component\":\"$COMPONENT\" } }"
curl --verbose -H "Accept: application/vnd.github.everest-preview+json" \
    -H "Authorization: token ${TOKEN}" \
    -H "Content-Type: application/json" \
    --request POST \
    --data "{\"event_type\": \"$EVENT_TYPE\", \"client_payload\": { \"dist\": \"$DIST\", \"version\": \"$version\", \"component\":\"$COMPONENT\" } }" \
    https://api.github.com/repos/${REPO}/dispatches


