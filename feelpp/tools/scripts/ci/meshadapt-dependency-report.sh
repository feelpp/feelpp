#!/usr/bin/env bash

set -u

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repository_root=${FEELPP_SOURCE_DIR:-}
if [[ -z "${repository_root}" ]]; then
    repository_root=$(git -C "${script_dir}" rev-parse --show-toplevel 2>/dev/null || true)
fi

print_first_line()
{
    local key=$1
    shift
    if command -v "$1" >/dev/null 2>&1; then
        printf '%s=%s\n' "${key}" "$("$@" 2>&1 | head -n 1)"
    else
        printf '%s=unavailable\n' "${key}"
    fi
}

print_package()
{
    local package=$1
    if command -v dpkg-query >/dev/null 2>&1; then
        local version
        version=$(dpkg-query -W -f='${Version}' "${package}" 2>/dev/null || true)
        printf 'package.%s=%s\n' "${package}" "${version:-not-installed}"
    elif command -v rpm >/dev/null 2>&1; then
        local version
        version=$(rpm -q --qf '%{VERSION}-%{RELEASE}' "${package}" 2>/dev/null || true)
        printf 'package.%s=%s\n' "${package}" "${version:-not-installed}"
    else
        printf 'package.%s=package-manager-unavailable\n' "${package}"
    fi
}

print_submodule()
{
    local name=$1
    local path=$2
    if [[ -d "${path}/.git" || -f "${path}/.git" ]]; then
        local commit describe
        commit=$(git -C "${path}" rev-parse HEAD 2>/dev/null || true)
        describe=$(git -C "${path}" describe --tags --always --dirty 2>/dev/null || true)
        printf '%s.commit=%s\n' "${name}" "${commit:-worktree-gitdir-unavailable}"
        printf '%s.describe=%s\n' "${name}" "${describe:-worktree-gitdir-unavailable}"
    else
        printf '%s=unavailable\n' "${name}"
    fi
}

if [[ -r /etc/os-release ]]; then
    os_pretty_name=$(sed -n 's/^PRETTY_NAME=//p' /etc/os-release | tr -d '"')
    printf 'os=%s\n' "${os_pretty_name}"
fi
printf 'kernel=%s\n' "$(uname -srmo)"
printf 'architecture=%s\n' "$(uname -m)"

print_first_line cc "${CC:-cc}" --version
print_first_line cxx "${CXX:-c++}" --version
print_first_line cmake cmake --version
print_first_line ninja ninja --version
if command -v ompi_info >/dev/null 2>&1; then
    print_first_line mpi ompi_info --version
else
    print_first_line mpi mpiexec --version
fi

for package in clang libscotch-dev libopenmpi-dev openmpi-bin cmake ninja-build; do
    print_package "${package}"
done

scotch_library=
if command -v "${CC:-cc}" >/dev/null 2>&1; then
    scotch_library=$("${CC:-cc}" -print-file-name=libscotch.so 2>/dev/null || true)
    if [[ "${scotch_library}" == "libscotch.so" || ! -r "${scotch_library}" ]]; then
        scotch_library=
    fi
fi
if command -v ldconfig >/dev/null 2>&1; then
    if [[ -z "${scotch_library}" ]]; then
        scotch_library=$(ldconfig -p 2>/dev/null | awk '$1 ~ /^libscotch-[0-9].*\.so$/ { print $NF; exit }')
    fi
fi
printf 'scotch.library=%s\n' "${scotch_library:-not-found}"
if [[ -n "${scotch_library}" ]] && command -v nm >/dev/null 2>&1; then
    if nm -D --defined-only "${scotch_library}" 2>/dev/null | awk '{print $3}' | grep -Fxq '_SCOTCHintSort2asc1'; then
        printf 'scotch.private_sort_symbol=exported\n'
    else
        printf 'scotch.private_sort_symbol=not-exported\n'
    fi
else
    printf 'scotch.private_sort_symbol=unknown\n'
fi

if [[ -n "${repository_root}" && -d "${repository_root}" ]]; then
    feelpp_commit=$(git -C "${repository_root}" rev-parse HEAD 2>/dev/null || true)
    printf 'feelpp.commit=%s\n' "${feelpp_commit:-worktree-gitdir-unavailable}"
    print_submodule mmg "${repository_root}/feelpp/contrib/mmg"
    print_submodule parmmg "${repository_root}/feelpp/contrib/parmmg"
    if grep -R -n -m 1 '_SCOTCHintSort2asc1' \
        "${repository_root}/feelpp/contrib/mmg/src/common" >/dev/null 2>&1; then
        printf 'mmg.private_sort_call=present\n'
    else
        printf 'mmg.private_sort_call=absent\n'
    fi
fi
