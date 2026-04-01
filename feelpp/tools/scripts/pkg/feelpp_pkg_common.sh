#! /bin/bash

set -eo pipefail

BUILDKITE_AGENT_NAME=${BUILDKITE_AGENT_NAME:-default}
BUILDKITE_BRANCH=${BUILDKITE_BRANCH:-${GITHUB_REF_NAME:-develop}}
BRANCH=${BRANCH:-${BUILDKITE_BRANCH}}

# default values
CHANNEL=${CHANNEL:-latest}
if [ "$BUILDKITE_BRANCH" = "develop" -o  "$BRANCH" = "develop" ]; then
    CHANNEL=latest
fi
if [ "$BUILDKITE_BRANCH" = "main" -o "$BRANCH" = "main" -o "$BUILDKITE_BRANCH" = "master" -o  "$BRANCH" = "master" ]; then
    CHANNEL=stable
fi 
DIST=${DIST:-noble}
case "$DIST" in
    focal|jammy|lunar|mantic|noble)
        FLAVOR=ubuntu
        ;;
    bullseye|bookworm|trixie|testing|sid)
        FLAVOR=debian
        ;;
    fedora-42)
        FLAVOR=fedora
        ;;
    *)
        echo "Unsupported DIST: $DIST" >&2
        exit 1
        ;;
esac



COMPONENT=${COMPONENT:-feelpp}
FEELPP_PKG_SCRIPT_DIR=${FEELPP_PKG_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
FEELPP_REPO_ROOT=${FEELPP_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}
FEELPP_PKG_DIR=${FEELPP_PKG_DIR:-${FEELPP_REPO_ROOT}/packaging/debian}
FEELPP_PBUILDER_DIR=${FEELPP_PBUILDER_DIR:-${FEELPP_REPO_ROOT}/packaging/pbuilder}
FEELPP_PBUILDER_CONFIG=${FEELPP_PBUILDER_CONFIG:-${FEELPP_PBUILDER_DIR}/pbuilderrc}
FEELPP_PBUILDER_SOURCE_HOOKDIR=${FEELPP_PBUILDER_SOURCE_HOOKDIR:-${FEELPP_PBUILDER_DIR}/hooks}
FEELPP_PBUILDER_AUTH=${FEELPP_PBUILDER_AUTH:-${FEELPP_PKG_SCRIPT_DIR}/feelpp_pkg_sudo_auth.sh}
FEELPP_PBUILDER_ROOT=${FEELPP_PBUILDER_ROOT:-${HOME}/pbuilder/chroots/${CHANNEL}}
FEELPP_PKG_JOB_ID=${FEELPP_PKG_JOB_ID:-${GITHUB_RUN_ID:-local}-${GITHUB_RUN_ATTEMPT:-0}-${GITHUB_JOB:-pkg}}
FEELPP_PKG_JOB_ROOT=${FEELPP_PKG_JOB_ROOT:-${HOME}/pbuilder/jobs/${FEELPP_PKG_JOB_ID}/${FLAVOR}/${DIST}}
FEELPP_PKG_LOCAL_REPO_DIR=${FEELPP_PKG_LOCAL_REPO_DIR:-${FEELPP_PKG_JOB_ROOT}/local-repo}
FEELPP_PKG_ARTIFACTS_DIR=${FEELPP_PKG_ARTIFACTS_DIR:-${FEELPP_PKG_JOB_ROOT}/artifacts}
FEELPP_PKG_RESULTS_DIR=${FEELPP_PKG_RESULTS_DIR:-${FEELPP_PKG_JOB_ROOT}/results}
FEELPP_PBUILDER_KEYRINGS_DIR=${FEELPP_PBUILDER_KEYRINGS_DIR:-${FEELPP_PKG_JOB_ROOT}/pbuilder/keyrings}
FEELPP_PBUILDER_HOOKDIR=${FEELPP_PBUILDER_HOOKDIR:-${FEELPP_PKG_JOB_ROOT}/pbuilder/hooks}
FEELPP_PBUILDER_MIRRORSITE=${FEELPP_PBUILDER_MIRRORSITE:-}
FEELPP_PBUILDER_OTHERMIRROR=${FEELPP_PBUILDER_OTHERMIRROR:-}
FEELPP_PKG_COMPONENT_DIR=${FEELPP_PKG_COMPONENT_DIR:-}
FEELPP_PKG_DIST_DIR=${FEELPP_PKG_DIST_DIR:-}

feelpp_pbuilder_root() {
    printf '%s\n' "${FEELPP_PBUILDER_ROOT}"
}

feelpp_pbuilder_basetgz() {
    local dist=${1:-$DIST}
    printf '%s/%s-base.tgz\n' "$(feelpp_pbuilder_root)" "$dist"
}

feelpp_pbuilder_basetgz_is_valid() {
    local base_tgz=${1:?base tarball path is required}

    if [ ! -f "${base_tgz}" ]; then
        return 1
    fi

    tar -tzf "${base_tgz}" >/dev/null 2>&1
}

prepare_feelpp_pkg_job_workspace() {
    mkdir -p "${FEELPP_PKG_LOCAL_REPO_DIR}" "${FEELPP_PKG_ARTIFACTS_DIR}" "${FEELPP_PKG_RESULTS_DIR}" "${FEELPP_PBUILDER_KEYRINGS_DIR}" "${FEELPP_PBUILDER_HOOKDIR}"
}

ensure_feelpp_pbuilder_runtime_assets() {
    if [ ! -d "${FEELPP_PBUILDER_HOOKDIR}" ] || [ -z "$(find "${FEELPP_PBUILDER_HOOKDIR}" -maxdepth 1 -type f -print -quit 2>/dev/null)" ]; then
        echo "pbuilder runtime hook directory is missing or empty: ${FEELPP_PBUILDER_HOOKDIR}" >&2
        echo "Run feelpp-pkg pbuilder prepare or feelpp-pkg build ... first." >&2
        exit 1
    fi

    if [ ! -d "${FEELPP_PBUILDER_KEYRINGS_DIR}" ] || [ -z "$(find "${FEELPP_PBUILDER_KEYRINGS_DIR}" -maxdepth 1 -type f -name '*.gpg' -print -quit 2>/dev/null)" ]; then
        echo "pbuilder runtime keyring directory is missing or empty: ${FEELPP_PBUILDER_KEYRINGS_DIR}" >&2
        echo "Run feelpp-pkg pbuilder prepare or feelpp-pkg build ... first." >&2
        exit 1
    fi
}

feelpp-pbuilder-dist() {
    local dist=$1
    local operation
    shift

    if [ $# -lt 1 ]; then
        echo "pbuilder operation is required" >&2
        exit 1
    fi

    if [ ! -f "${FEELPP_PBUILDER_CONFIG}" ]; then
        echo "pbuilder config not found: ${FEELPP_PBUILDER_CONFIG}" >&2
        exit 1
    fi

    if [ -z "${FEELPP_PBUILDER_MIRRORSITE}" ]; then
        echo "pbuilder mirrorsite is not configured" >&2
        echo "Run feelpp-pkg build ... or set FEELPP_PBUILDER_MIRRORSITE." >&2
        exit 1
    fi

    if [ -z "${FEELPP_PBUILDER_OTHERMIRROR}" ]; then
        echo "pbuilder othermirror set is not configured" >&2
        echo "Run feelpp-pkg build ... or set FEELPP_PBUILDER_OTHERMIRROR." >&2
        exit 1
    fi

    mkdir -p "$(feelpp_pbuilder_root)"
    ensure_feelpp_pbuilder_runtime_assets

    operation=$1
    shift

    env \
        PBUILDFOLDER="$(feelpp_pbuilder_root)" \
        PBUILDAUTH="${FEELPP_PBUILDER_AUTH}" \
        FEELPP_PBUILDER_KEYRINGS_DIR="${FEELPP_PBUILDER_KEYRINGS_DIR}" \
        FEELPP_PBUILDER_BINDMOUNTS="${FEELPP_PKG_LOCAL_REPO_DIR}" \
        MIRRORSITE="${FEELPP_PBUILDER_MIRRORSITE}" \
        OTHERMIRROR="${FEELPP_PBUILDER_OTHERMIRROR}" \
        pbuilder-dist "$dist" "$operation" --updates-only \
        --configfile "${FEELPP_PBUILDER_CONFIG}" \
        --hookdir "${FEELPP_PBUILDER_HOOKDIR}" \
        "$@"
}

prepare_feelpp_pbuilder_base() {
    local dist=${1:-$DIST}
    local base_tgz

    base_tgz=$(feelpp_pbuilder_basetgz "$dist")
    if [ -f "${base_tgz}" ] && ! feelpp_pbuilder_basetgz_is_valid "${base_tgz}"; then
        echo "--- removing invalid pbuilder base at ${base_tgz}"
        rm -f "${base_tgz}"
    fi

    if [ ! -f "${base_tgz}" ]; then
        echo "--- creating pbuilder base at ${base_tgz}"
        feelpp-pbuilder-dist "$dist" create
    else
        echo "--- updating pbuilder base at ${base_tgz}"
        feelpp-pbuilder-dist "$dist" update
    fi
}

prepare_feelpp_packaging_tree() {
    if [ -z "${FEELPP_PKG_COMPONENT_DIR}" ]; then
        FEELPP_PKG_COMPONENT_DIR="${FEELPP_PKG_DIR}/${COMPONENT}"
    fi
    if [ -z "${FEELPP_PKG_DIST_DIR}" ]; then
        FEELPP_PKG_DIST_DIR="${FEELPP_PKG_COMPONENT_DIR}/${DIST}"
    fi

    if [ ! -d "${FEELPP_PKG_COMPONENT_DIR}" ]; then
        echo "Packaging metadata not found for component ${COMPONENT}: ${FEELPP_PKG_COMPONENT_DIR}" >&2
        exit 1
    fi

    if [ ! -d "${FEELPP_PKG_DIST_DIR}" ]; then
        echo "Packaging metadata not found for ${COMPONENT}/${DIST}: ${FEELPP_PKG_DIST_DIR}" >&2
        echo "Available distros for ${COMPONENT}:" >&2
        ls -1 "${FEELPP_PKG_COMPONENT_DIR}" >&2
        exit 1
    fi

    echo "--- using in-tree packaging metadata from ${FEELPP_PKG_DIST_DIR}"
}
