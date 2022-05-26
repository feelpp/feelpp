 #!/bin/bash

set -euo pipefail

FEELPP_DIR=${FEELPP_DIR:-$PWD}

source $(dirname $0)/common.sh

## Prints supported image combinations in the following format
## {image_name} {base_image} {distro} {feelpp_version} {docker_version} {extra tags...}

DEBIAN_VERSIONS=(
    10
    testing
    sid
)

LATEST_DEBIAN=${DEBIAN_VERSIONS[${#DEBIAN_VERSIONS[@]} - 2]}

UBUNTU_VERSIONS=(
    22.04
    20.04
    19.10
    19.04
    18.04
    #16.04
    17.04
    16.10
    16.04

)

LATEST_UBUNTU=${UBUNTU_VERSIONS[${#UBUNTU_VERSIONS[@]} - 1]}

DISTROS=(
  debian
  ubuntu
)

FEELPP_BRANCHES=(
    #master
    #develop
    $1
    #${BRANCHTAG}
)

# Returns the major version for a given X.X.X docker version
docker_major_version() {
    #cut -d. -f1-2 <<< $1
    cut -d. -f1 <<< $1
}

image_name() {
  local release="$1"
  local distro="$2"
  printf "%s-%s" "$release" "$distro" | sed -e 's/^stable-//g'
}

version=${2:-${FEELPP_VERSION}}

for feelpp_branch in ${FEELPP_BRANCHES[*]} ; do
  distro="debian"
  docker="n/a"
  feelpp_branch_version="${feelpp_branch}-${version}"
  feelpp_version="${feelpp_branch_version}"
  for os_version in ${DEBIAN_VERSIONS[*]} ; do
    printf "%s-%s" \
      $(image_name "$feelpp_branch_version" "$distro") $os_version 
#      "n/a" \
#      $feelpp_branch_version $feelpp_version $distro $os_version  \
#      $(image_name "$feelpp_branch_version" "$distro") $(docker_major_version $os_version)

#    if [[ $os_version == $LATEST_DEBIAN ]] ; then
#      printf " %s" $(image_name "$feelpp_version" "$distro")

      #printf " %s" $(image_name "$feelpp_version" "$distro")
#    fi

    echo #newline
  done
done

for feelpp_branch in ${FEELPP_BRANCHES[*]} ; do
  distro="ubuntu"
  docker="n/a"
  feelpp_branch_version="${feelpp_branch}-${version}"
  feelpp_version="${feelpp_branch}"
#  printf "%s %s %s %s %s\n" $(image_name "$feelpp_version" "$distro") "n/a" "$distro" "$feelpp_version" "$docker"
  for os_version in ${UBUNTU_VERSIONS[*]} ; do
      #    printf "%s-%s %s %s %s %s %s %s-%s" \
      printf "%s-%s" \
      $(image_name "$feelpp_branch_version" "$distro") $os_version 
#      $(image_name "$feelpp_branch_version" "$distro") \
#      $feelpp_branch_version $feelpp_version $distro $os_version  \
#      $(image_name "$feelpp_branch_version" "$distro") $(docker_major_version $os_version)

    if [[ $os_version == $LATEST_UBUNTU ]] ; then
        # We also want to give the ubuntu distro the official
        # feelpp/<container>:[latest,edge,beta] tags

        printf " %s" $(sed -e 's/develop/latest/g' <<< $feelpp_version)
        printf " %s" $(sed -e 's/master/stable/g' <<< $feelpp_version)

        # tag with version only with master and develop
        # master should take over develop once a release is done
        if [ ${feelpp_branch} = "develop" -o ${feelpp_branch} = "master" ]; then
            printf " %s" "${version}"
        fi

#        printf " %s" $(image_name "$feelpp_version" "$distro")
    fi

    echo #newline
  done
done
