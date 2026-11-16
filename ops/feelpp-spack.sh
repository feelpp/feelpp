# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later
# Source once from Bash to enable environments in the current shell.

_feelpp_ops_file="${BASH_SOURCE[0]}"
case "$_feelpp_ops_file" in
    /*) ;;
    *) _feelpp_ops_file="$PWD/$_feelpp_ops_file" ;;
esac
_feelpp_ops_directory="${_feelpp_ops_file%/*}"
_feelpp_ops_source="$_feelpp_ops_directory/src"
export FEELPP_REPO_ROOT="${_feelpp_ops_directory%/*}"
if [ -n "${FEELPP_OPS_PYTHON:-}" ]; then
    _feelpp_ops_python="$FEELPP_OPS_PYTHON"
elif [ -x "$FEELPP_REPO_ROOT/ops/.venv/bin/python" ]; then
    _feelpp_ops_python="$FEELPP_REPO_ROOT/ops/.venv/bin/python"
elif [ -x "$FEELPP_REPO_ROOT/.venv-ops/bin/python" ]; then
    _feelpp_ops_python="$FEELPP_REPO_ROOT/.venv-ops/bin/python"
else
    _feelpp_ops_python="$(command -v python3)"
fi

# Use this checkout's CLI even when its console entry point is not installed.
_fpp_spack_cli()
{
    PYTHONPATH="${_feelpp_ops_source}${PYTHONPATH:+:$PYTHONPATH}" \
        "$_feelpp_ops_python" -m feelpp.pkg spack "$@"
}

# Only enable/disable evaluate shell changes. Activate retains its script output.
fpp-spack()
{
    local _feelpp_spack_script
    case "${1:-}" in
        enable)
            shift
            _feelpp_spack_script="$(_fpp_spack_cli activate "$@")" || return
            eval "$_feelpp_spack_script"
            ;;
        disable)
            if [ "$#" -ne 1 ]; then
                printf '%s\n' 'usage: fpp-spack disable' >&2
                return 2
            fi
            if [ -n "${SPACK_ENV:-}" ]; then
                _feelpp_spack_script="$(command "$SPACK_ROOT/bin/spack" env deactivate --sh)" || return
                eval "$_feelpp_spack_script"
            fi
            ;;
        *) _fpp_spack_cli "$@" ;;
    esac
}

# Read whole candidate lines so environment names containing spaces stay intact.
# The pinned ops interpreter keeps completion working after Spack changes PATH.
_fpp_spack_complete()
{
    local _feelpp_candidate
    COMPREPLY=()
    while IFS= read -r _feelpp_candidate; do
        COMPREPLY+=("$_feelpp_candidate")
    done < <(
        PYTHONPATH="${_feelpp_ops_source}${PYTHONPATH:+:$PYTHONPATH}" \
            "$_feelpp_ops_python" -m feelpp.pkg.backends.spack.completion \
            "${COMP_WORDS[@]:1:COMP_CWORD}" 2>/dev/null
    )
}
complete -o filenames -F _fpp_spack_complete fpp-spack
