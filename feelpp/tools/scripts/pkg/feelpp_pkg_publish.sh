#!/bin/bash

set -eo pipefail
if [ "${FEELPP_PKG_DEBUG_SH:-false}" = "true" ]; then
    set -x
fi
source $(dirname $0)/feelpp_pkg_common.sh

repo_name="feelpp-${DIST}-${CHANNEL}"
snapshot_id=${FEELPP_PKG_SNAPSHOT_ID:-${GITHUB_RUN_ID:-local}-${GITHUB_RUN_ATTEMPT:-0}-$(date -u +%Y%m%d%H%M%S)}
snapshot_name="${repo_name}-snapshot-${snapshot_id}"
publish_distribution=${FEELPP_APTLY_PUBLISH_DISTRIBUTION:-${DIST}}
publish_component=${FEELPP_APTLY_PUBLISH_COMPONENT:-${CHANNEL}}
publish_target=${FEELPP_APTLY_PUBLISH_TARGET:-s3:apt.feelpp.org:${FLAVOR}/${DIST}}
passphrase_file=""

cleanup() {
    if [ -n "${passphrase_file}" ] && [ -f "${passphrase_file}" ]; then
        rm -f "${passphrase_file}"
    fi
}

trap cleanup EXIT

aptly_args=()
publish_args=(-force-overwrite)

if [ -n "${FEELPP_APTLY_CONFIG:-}" ]; then
    aptly_args+=(-config="${FEELPP_APTLY_CONFIG}")
fi

if [ "${FEELPP_APTLY_SKIP_SIGNING:-false}" = "true" ]; then
    publish_args+=(-skip-signing)
else
    if [ -n "${GPG_PASSPHRASE:-}" ]; then
        passphrase_file=$(mktemp "${FEELPP_PKG_JOB_ROOT}/aptly-passphrase.XXXXXX")
        chmod 600 "${passphrase_file}"
        printf '%s' "${GPG_PASSPHRASE}" > "${passphrase_file}"
        publish_args+=(-batch -passphrase-file="${passphrase_file}")
    fi
    if [ -n "${GPG_KEY:-}" ]; then
        publish_args+=(-gpg-key="${GPG_KEY}")
    fi
fi

run_aptly() {
    aptly "${aptly_args[@]}" "$@"
}

if ! run_aptly repo show "${repo_name}" >/dev/null 2>&1; then
    run_aptly repo create -distribution="${publish_distribution}" -component="${publish_component}" "${repo_name}"
fi

if [ -n "${FEELPP_PKG_PUBLISH_INPUT_DIR:-}" ]; then
    if [ ! -d "${FEELPP_PKG_PUBLISH_INPUT_DIR}" ]; then
        echo "Publish input directory not found: ${FEELPP_PKG_PUBLISH_INPUT_DIR}" >&2
        exit 1
    fi

    if find "${FEELPP_PKG_PUBLISH_INPUT_DIR}" -type f \( -name '*.deb' -o -name '*.udeb' \) | grep -q .; then
        run_aptly repo add -force-replace "${repo_name}" "${FEELPP_PKG_PUBLISH_INPUT_DIR}"
    else
        echo "No binary packages found in ${FEELPP_PKG_PUBLISH_INPUT_DIR}, skipping aptly repo add"
    fi
fi

run_aptly snapshot create "${snapshot_name}" from repo "${repo_name}"

if run_aptly publish show "${publish_distribution}" "${publish_target}" >/dev/null 2>&1; then
    run_aptly publish switch "${publish_args[@]}" -component="${publish_component}" "${publish_distribution}" "${publish_target}" "${snapshot_name}"
else
    run_aptly publish snapshot "${publish_args[@]}" -distribution="${publish_distribution}" -component="${publish_component}" "${snapshot_name}" "${publish_target}"
fi

echo "Published snapshot: ${snapshot_name}"
