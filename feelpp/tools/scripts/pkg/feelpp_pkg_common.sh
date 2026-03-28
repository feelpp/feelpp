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
FEELPP_PKG_DIR=${FEELPP_PKG_DIR:-feelpp.pkg}
FEELPP_PKG_REPO=${FEELPP_PKG_REPO:-https://github.com/feelpp/feelpp.pkg.git}
FEELPP_PKG_REF=${FEELPP_PKG_REF:-}

# Define the function
feelpp-pbuilder-dist() {
    local dist=$1
    shift  # Shift arguments to pass any additional parameters to pbuilder-dist

    if [ "$dist" = "toto" ]; then
        # Handle the 'noble' distribution with an automatic 'y' response
        {
            echo y
            cat
        } | pbuilder-dist "$dist" "$@"
    else
        # Handle other distributions normally
        pbuilder-dist "$dist" "$@"
    fi
}

prepare_feelpp_pkg_checkout() {
    local cloned=0

    if [ ! -d "${FEELPP_PKG_DIR}/.git" ]; then
        echo "--- cloning feelpp.pkg from ${FEELPP_PKG_REPO}"
        git clone -q "${FEELPP_PKG_REPO}" "${FEELPP_PKG_DIR}"
        cloned=1
    else
        echo "--- using existing ${FEELPP_PKG_DIR} checkout"
    fi

    if [ -n "${FEELPP_PKG_REF}" ]; then
        echo "--- pinning feelpp.pkg to ${FEELPP_PKG_REF}"
        (
            cd "${FEELPP_PKG_DIR}"
            git fetch -q --all --tags
            git checkout -q "${FEELPP_PKG_REF}"
        )
    elif [ "$cloned" -eq 1 ]; then
        echo "--- FEELPP_PKG_REF not set; using repository default branch"
    fi

    echo "--- feelpp.pkg revision: $(git -C "${FEELPP_PKG_DIR}" rev-parse HEAD)"
}
