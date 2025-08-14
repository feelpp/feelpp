#!/usr/bin/env bash
set -euo pipefail
# Downstream may close early (head/awk exit); don't treat EPIPE as fatal.
trap '' PIPE

# --- Config / Inputs ---------------------------------------------------------
FEELPP_DIR=${FEELPP_DIR:-$PWD}
source "$(dirname "$0")/common.sh"

# Usage: script <feelpp_branch> [feelpp_version]
if [[ $# -lt 1 ]]; then
  echo "usage: $(basename "$0") <feelpp_branch> [feelpp_version]" >&2
  exit 1
fi

FEELPP_BRANCH="$1"
FEELPP_VERSION_INPUT="${2:-${FEELPP_VERSION:-}}"
if [[ -z "${FEELPP_VERSION_INPUT}" ]]; then
  echo "FEELPP_VERSION not provided (arg #2 or env FEELPP_VERSION)." >&2
  exit 1
fi

# Versions per distro
DEBIAN_VERSIONS=(13 12 11 testing sid)
UBUNTU_VERSIONS=(24.04 23.10 22.04 20.04)
FEDORA_VERSIONS=(42)

# --- Helpers -----------------------------------------------------------------
docker_major_version() {
  # If you ever need the major-only again:
  # cut -d. -f1 <<< "$1"
  cut -d. -f1 <<< "$1"
}

image_name() {
  # {release}-{distro} with 'stable-' prefix stripped
  local release="$1" distro="$2"
  printf "%s-%s" "$release" "$distro" | sed -e 's/^stable-//g'
}

latest_of() {
  # "Latest" = last element of an array passed by name
  local -n _arr="$1"
  echo "${_arr[${#_arr[@]}-1]}"
}

print_debian_lines() {
  local branch="$1" version="$2"
  local distro="debian"
  local branch_version="${branch}-${version}"

  for os_version in "${DEBIAN_VERSIONS[@]}"; do
    printf "%s-%s\n" "$(image_name "$branch_version" "$distro")" "$os_version"
  done
}

print_ubuntu_lines() {
  local branch="$1" version="$2"
  local distro="ubuntu"
  local branch_version="${branch}-${version}"
  local latest_ubuntu
  latest_ubuntu="$(latest_of UBUNTU_VERSIONS)"

  for os_version in "${UBUNTU_VERSIONS[@]}"; do
    # Base tag
    line="$(printf "%s-%s" "$(image_name "$branch_version" "$distro")" "$os_version")"

    # On the latest Ubuntu only, append extra tags:
    if [[ "$os_version" == "$latest_ubuntu" ]]; then
      # latest/stable mapping derived from branch
      # develop -> latest, master -> stable
      line+=" $(sed -e 's/develop/latest/g' <<< "$branch")"
      line+=" $(sed -e 's/master/stable/g'  <<< "$branch")"

      # If branch is develop or master, also tag with numeric FEELPP version
      if [[ "$branch" == "develop" || "$branch" == "master" ]]; then
        line+=" ${version}"
      fi
    fi

    printf "%s\n" "$line"
  done
}

# If you want Fedora printed similarly later, this is ready.
print_fedora_lines() {
  local branch="$1" version="$2"
  local distro="fedora"
  local branch_version="${branch}-${version}"

  for os_version in "${FEDORA_VERSIONS[@]}"; do
    printf "%s-%s\n" "$(image_name "$branch_version" "$distro")" "$os_version"
  done
}

# --- Output ------------------------------------------------------------------
# Debian block (no special tag expansions)
print_debian_lines "$FEELPP_BRANCH" "$FEELPP_VERSION_INPUT"

# Ubuntu block (adds extra tags on latest)
print_ubuntu_lines "$FEELPP_BRANCH" "$FEELPP_VERSION_INPUT"

# If/when you want Fedora lines too, uncomment:
print_fedora_lines "$FEELPP_BRANCH" "$FEELPP_VERSION_INPUT"