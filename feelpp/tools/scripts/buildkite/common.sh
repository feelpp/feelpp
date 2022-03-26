#! /bin/bash

set -eo pipefail

FEELPP_DIR=${FEELPP_DIR:-/usr/local}
if [ -d ${FEELPP_DIR}/share/feelpp/feel/cmake/modules/ ]; then 
    FEELPP_CMAKE_DIR=${FEELPP_DIR}/share/feelpp/feel/cmake/modules;
else
    FEELPP_CMAKE_DIR=. #feelpp/cmake/modules
fi
if [ -d ${FEELPP_DIR}/share/feelpp/scripts/ ]; then 
    FEELPP_SCRIPTS_DIR=${FEELPP_DIR}/share/feelpp/scripts;
else
    FEELPP_SCRIPTS_DIR=feelpp/tools/scripts/buildkite
fi

function get_field(){
    local fname="$1"
    local field="$2"
    ret=$(cat $fname | grep $field | awk 'BEGIN{FS=" "}{print $2}' | awk 'BEGIN{FS="\""}{print $2}')
    printf "%s" "${ret}"
}
function get_version() {
    vfile=${FEELPP_CMAKE_DIR}/feelpp.version.cmake
    IFS=''
    major=$(get_field "$vfile" "VERSION_MAJOR ")
    minor=$(get_field "$vfile" "VERSION_MINOR ")
    micro=$(get_field "$vfile" "VERSION_MICRO ")
    prerelease=$(get_field "$vfile" "VERSION_PRERELEASE ")
    printf "v%s" "${major}.${minor}.${micro}${prerelease}"
}

tag_from_target() {
    splitfrom=(`echo "$1" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}

    ${FEELPP_SCRIPTS_DIR}/list.sh $2 $3 | grep "$2-$3-${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        image=${tokens[0]}
        printf "%s" "${image}"
    done
}
tag_from_os() {
    splitfrom=(`echo "$1" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}
    printf "%s" "${fromos}-${fromtag}"
    #    ${FEELPP_SCRIPTS_DIR}/list.sh $2 $3 | grep "${fromos}-${fromtag}"  | while read line ; do
    #    tokens=($line)
    #    image=${tokens[0]}
    #    printf "%s" "${image}"
    #done
}

extratags_from_target() {
    splitfrom=(`echo "$1" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}
    
    ${FEELPP_SCRIPTS_DIR}/list.sh $2 $3 | grep "$2-$3-${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        extratags=${tokens[@]:1:3}
        printf "%s" "${extratags}" 
    done
}

# Combines a dockerfile template with a generated FROM line
dockerfile_from() {
    local dockerfile="$1"
    local from="$2"
    printf 'FROM %s\n%s' "$from" "$(<$dockerfile)"
}

FEELPP_VERSION=$(get_version)
#echo "Feel++ version: $FEELPP_VERSION"
