#!/usr/bin/env bash

set -euo pipefail

env_args=()
while [ $# -gt 0 ]; do
    case "$1" in
        *=*)
            env_args+=("$1")
            shift
            ;;
        *)
            break
            ;;
    esac
done

if [ $# -eq 0 ]; then
    echo "feelpp_pkg_sudo_auth.sh: missing command" >&2
    exit 1
fi

if [ "$(id -u)" -eq 0 ]; then
    exec env "${env_args[@]}" "$@"
fi

exec sudo env "${env_args[@]}" "$@"
