#! /bin/bash

set -eo pipefail

function get_field(){
    local fname="$1"
    local field="$2"
    ret=$(cat $fname | grep $field | awk 'BEGIN{FS=" "}{print $2}' | awk 'BEGIN{FS="\""}{print $2}')
    printf "%s" "${ret}"
}
function get_version() {
    IFS=''
    major=$(get_field "cmake/modules/feelpp.version.cmake" "VERSION_MAJOR ")
    minor=$(get_field "cmake/modules/feelpp.version.cmake" "VERSION_MINOR ")
    micro=$(get_field "cmake/modules/feelpp.version.cmake" "VERSION_MICRO ")
    prerelease=$(get_field "cmake/modules/feelpp.version.cmake" "VERSION_PRERELEASE ")
    printf "v%s" "${major}.${minor}.${micro}${prerelease}"
}

tag_from_target() {
    splitfrom=(`echo "$TARGET" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}

    tools/scripts/buildkite/list.sh | grep "${BRANCHTAG}-${version}-${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        image=${tokens[0]}
        printf "%s" "${image}"
    done
}
extratags_from_target() {
    splitfrom=(`echo "$TARGET" | tr ":" "\n"`)
    fromos=${splitfrom[0]}
    fromtag=${splitfrom[1]}
    
    tools/scripts/buildkite/list.sh | grep "${BRANCHTAG}-${version}-${fromos}-${fromtag}"  | while read line ; do
        tokens=($line)
        extratags=${tokens[@]:5}
        printf "%s" "${extratags}" 
    done
}

# Combines a dockerfile template with a generated FROM line
dockerfile_from() {
    local dockerfile="$1"
    local from="$2"
    printf 'FROM %s\n%s' "$from" "$(<$dockerfile)"
}

