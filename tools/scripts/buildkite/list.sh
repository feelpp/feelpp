 #!/bin/bash

set -euo pipefail

source $(dirname $0)/common.sh

## Prints supported image combinations in the following format
## {image_name} {base_image} {distro} {feelpp_version} {docker_version} {extra tags...}

DEBIAN_VERSIONS=(
    testing
    sid
)

LATEST_DEBIAN=${DEBIAN_VERSIONS[${#DEBIAN_VERSIONS[@]} - 2]}

UBUNTU_VERSIONS=(
    #16.04
    17.04
    16.10
)

LATEST_UBUNTU=${UBUNTU_VERSIONS[${#UBUNTU_VERSIONS[@]} - 1]}

DISTROS=(
  debian
  ubuntu
)

FEELPP_BRANCHES=(
    #master
    #develop
    ${BRANCHTAG}
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

version=$(get_version)

for feelpp_branch in ${FEELPP_BRANCHES[*]} ; do
  distro="debian"
  docker="n/a"
  feelpp_branch_version="${feelpp_branch}-${version}"
  feelpp_version="${feelpp_branch_version}"
  for os_version in ${DEBIAN_VERSIONS[*]} ; do
    printf "%s-%s %s %s %s %s %s %s-%s" \
      $(image_name "$feelpp_branch_version" "$distro") $os_version \
      "n/a" \
      $distro $feelpp_version $os_version $feelpp_branch_version \
      $(image_name "$feelpp_branch_version" "$distro") $(docker_major_version $os_version)

    if [[ $os_version == $LATEST_DEBIAN ]] ; then
      printf " %s" $(image_name "$feelpp_version" "$distro")

      #printf " %s" $(image_name "$feelpp_version" "$distro")
    fi

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
    printf "%s-%s %s %s %s %s %s %s-%s" \
      $(image_name "$feelpp_branch_version" "$distro") $os_version \
      $(image_name "$feelpp_branch_version" "$distro") \
      $distro $feelpp_version $os_version $feelpp_branch_version \
      $(image_name "$feelpp_branch_version" "$distro") $(docker_major_version $os_version)

    if [[ $os_version == $LATEST_UBUNTU ]] ; then
        # We also want to give the ubuntu distro the official
        # feelpp/<container>:[latest,edge,beta] tags

        printf " %s" $(sed -e 's/develop/latest/g' <<< $feelpp_version)
        printf " %s" $(sed -e 's/master/stable/g' <<< $feelpp_version)


        printf " %s" $(image_name "$feelpp_version" "$distro")
    fi

    echo #newline
  done
done
