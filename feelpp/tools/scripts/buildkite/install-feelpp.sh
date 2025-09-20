#!/usr/bin/env bash
# build-component.sh
#
# Usage:
#   ./build-component.sh [component]
# Env:
#   DRY_RUN=1           # show what would happen, do not execute
#   DEBUG=1             # enable bash -x tracing
#   TARGET              # e.g. "debian:13", required
#   BRANCH              # e.g. "develop" or "master", required
#   CC / CXX            # compiler tags
#   JOBS                # build jobs
#   CONFIGURE_FLAGS     # extra configure flags
#   CMAKE_FLAGS         # extra cmake flags
#   FEELPP_GITHUB_TOKEN # optional
#   FEELPP_GIRDER_API_KEY # optional

set -euo pipefail

# ---- debug / tracing ---------------------------------------------------------
DEBUG="${DEBUG:-0}"
if [[ "$DEBUG" == "1" ]]; then set -x; fi

# ---- dry-run helper ----------------------------------------------------------
DRY_RUN="${DRY_RUN:-0}"
run() {
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "[DRY-RUN] $*"
  else
    "$@"
  fi
}

say() { printf -- "%s\n" "$*"; }

# ---- script location & helpers ----------------------------------------------
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# We expect common.sh to define:
#   - dockerfile_from <template> <from>
#   - tag_from_target <TARGET> <BRANCHTAG> <FEELPP_VERSION>
#   - tag_from_os <TARGET> <BRANCHTAG> <FEELPP_VERSION>
#   - extratags_from_target <TARGET> <BRANCHTAG> <FEELPP_VERSION>
#   - get_version
# shellcheck source=/dev/null
source "${script_dir}/common.sh"

# ---- args --------------------------------------------------------------------
component="${1:-base}"

# ---- required env (give gentle defaults to pass --dry-run gracefully) -------
TARGET="${TARGET:-debian:13}"
BRANCH="${BRANCH:-develop}"
CXX="${CXX:-clang++}"
CC="${CC:-clang}"
JOBS="${JOBS:-$(nproc || echo 4)}"
CONFIGURE_FLAGS="${CONFIGURE_FLAGS:-}"
CMAKE_FLAGS="${CMAKE_FLAGS:-}"
FEELPP_GITHUB_TOKEN="${FEELPP_GITHUB_TOKEN:-}"
FEELPP_GIRDER_API_KEY="${FEELPP_GIRDER_API_KEY:-}"

# ---- compute tags ------------------------------------------------------------
say '--- clone/pull feelpp/docker'
if [[ -d "${script_dir}/docker" ]]; then
  run bash -lc "cd '${script_dir}/docker' && git pull --ff-only"
else
  run git clone --depth=1 https://github.com/feelpp/docker "${script_dir}/docker"
fi

BRANCHTAG="$(echo "${BRANCH}" | sed -e 's#/#-#g')"
FEELPP_VERSION="$(get_version)"
tag_compiler="$(echo "${CC}" | sed -e 's/-//g')"

base_tag="$(tag_from_target "${TARGET}" "${BRANCHTAG}" "${FEELPP_VERSION}")"
if [[ "${tag_compiler}" == gcc* ]]; then
  tag="${base_tag}-${tag_compiler}"
else
  tag="${base_tag}"
fi

tagos="$(tag_from_os "${TARGET}" "${BRANCHTAG}" "${FEELPP_VERSION}")"

# Optional debug suffix based on CI slug
BUILDKITE_PIPELINE_SLUG="${BUILDKITE_PIPELINE_SLUG:-}"
if [[ "${BUILDKITE_PIPELINE_SLUG}" == "feelpp-debug" ]]; then
  tag="${tag}-debug"
fi

# Image selection
image="feelpp-${component}"
<<<<<<< HEAD
if [ "${component}" = "feelpp" ] ; then
#    tag=$(tag_from_os $TARGET $BRANCHTAG $FEELPP_VERSION)
    image="feelpp"
elif [ "${component}" = "feelpp-core" ] ; then
    image="feelpp"    
fi
if [ "${component}" = "feelpp-python" ] ; then
#    tag=$(tag_from_os $TARGET $BRANCHTAG $FEELPP_VERSION)
    image="feelpp-python"
fi
echo "--- Building ${image}:${tag}"
=======
case "${component}" in
  feelpp)         image="feelpp" ;;
  feelpp-python|python)
                   image="feelpp-python" ;;
esac
>>>>>>> origin/develop

say "--- Building ${image}:${tag}"

<<<<<<< HEAD
if [ "${component}" = "feelpp" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-env:${tagos}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "feelpp-core" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-env:${tagos}" > docker/${image}/dockerfile.tmp    
elif [ "${component}" = "toolboxes" -o "${component}" = "testsuite" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp:${tag}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "mor" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-toolboxes:${tag}" > docker/${image}/dockerfile.tmp
elif [ "${component}" = "feelpp-python" -o "${component}" = "python" ] ; then
    dockerfile_from "docker/${image}/Dockerfile.template" "ghcr.io/feelpp/feelpp-mor:${tag}" > docker/${image}/dockerfile.tmp
=======
# ---- choose base image per component & generate dockerfile.tmp ---------------
tmp_df="${script_dir}/docker/${image}/dockerfile.tmp"
template="${script_dir}/docker/${image}/Dockerfile.template"

case "${component}" in
  feelpp)
    base_from="ghcr.io/feelpp/feelpp-env:${tagos}"
    ;;
  toolboxes|testsuite)
    base_from="ghcr.io/feelpp/feelpp:${tag}"
    ;;
  mor)
    base_from="ghcr.io/feelpp/feelpp-toolboxes:${tag}"
    ;;
  feelpp-python|python)
    base_from="ghcr.io/feelpp/feelpp-mor:${tag}"
    ;;
  *)
    base_from="ghcr.io/feelpp/feelpp-toolboxes:${tag}"
    ;;
esac

# Auto-detect platforms + description
ARCHES="$(arches_for_target "${TARGET}")"
DESCRIPTION="$(description_for "${image}" "${TARGET}")"

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[DRY-RUN] Generate ${tmp_df} from ${template} (BASE ${base_from})"
  dockerfile_from "${template}" "${base_from}" "${DESCRIPTION}"
>>>>>>> origin/develop
else
  mkdir -p "$(dirname "${tmp_df}")"
  dockerfile_from "${template}" "${base_from}" "${DESCRIPTION}" > "${tmp_df}"
fi

# ---- ctest flags per component ----------------------------------------------
case "${component}" in
  feelpp)          CTEST_FLAGS="-R feelpp_qs_ -T test --no-compress-output" ;;
  toolboxes)       CTEST_FLAGS="-R feelpp_toolbox_ -T test --no-compress-output --output-on-failure" ;;
  testsuite)       CTEST_FLAGS="-R feelpp_test_ -T test --no-compress-output --output-on-failure" ;;
  feelpp-python|python)
                   CTEST_FLAGS="-R feelpp -T test --no-compress-output --output-on-failure" ;;
  *)
                   CTEST_FLAGS="-T test --no-compress-output --output-on-failure" ;;
esac

# ---- docker build ------------------------------------------------------------


run docker build \
  --pull \
  --tag="ghcr.io/feelpp/${image}:${tag}" \
  --build-arg="FEELPP_GITHUB_TOKEN=${FEELPP_GITHUB_TOKEN}" \
  --build-arg="FEELPP_GIRDER_API_KEY=${FEELPP_GIRDER_API_KEY}" \
  --build-arg="BUILD_JOBS=${JOBS}" \
  --build-arg="BRANCH=${BRANCH}" \
  --build-arg="CXX=${CXX}" \
  --build-arg="CC=${CC}" \
  --build-arg="CONFIGURE_FLAGS=${CONFIGURE_FLAGS}" \
  --build-arg="CMAKE_FLAGS=${CMAKE_FLAGS}" \
  --build-arg="CTEST_FLAGS=${CTEST_FLAGS}" \
  --no-cache=true \
  -f "${tmp_df}" \
  "${script_dir}/docker/${image}"

# ---- extra tags --------------------------------------------------------------
say "--- Tagging ghcr.io/feelpp/${image}:${tag}"
read -r -a extra_tags <<< "$(extratags_from_target "${TARGET}" "${BRANCHTAG}" "${FEELPP_VERSION}")"

for tagalias in "${extra_tags[@]}"; do
  [[ -z "${tagalias}" ]] && continue
  say "Tagging ghcr.io/feelpp/${image}:${tag} as ghcr.io/feelpp/${image}:${tagalias}"
  run docker tag "ghcr.io/feelpp/${image}:${tag}" "ghcr.io/feelpp/${image}:${tagalias}"
done

# ---- release step ------------------------------------------------------------
if [[ "$DRY_RUN" == "1" ]]; then
  say "[DRY-RUN] Would call release.sh for ${image}"
else
  # shellcheck source=/dev/null
  source "${script_dir}/release.sh" -- "${image}"
fi

say "--- Done."