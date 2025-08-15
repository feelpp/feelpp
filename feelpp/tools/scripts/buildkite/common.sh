#!/usr/bin/env bash
set -euo pipefail

# ------------ Global switches ------------
DRY_RUN="${DRY_RUN:-0}"

# ------------ Config / defaults ------------
FEELPP_DIR="${FEELPP_DIR:-/usr/local}"
DEFAULT_TARGET="${DEFAULT_TARGET:-ubuntu:24.04}"

# ------------ Locate CMake + scripts dirs ------------
if [[ -d "${FEELPP_DIR}/share/feelpp/feel/cmake/modules" ]]; then
  FEELPP_CMAKE_DIR="${FEELPP_DIR}/share/feelpp/feel/cmake/modules"
else
  FEELPP_CMAKE_DIR="."
fi

if [[ -d "${FEELPP_DIR}/share/feelpp/scripts" ]]; then
  FEELPP_SCRIPTS_DIR="${FEELPP_DIR}/share/feelpp/scripts"
else
  FEELPP_SCRIPTS_DIR="feelpp/tools/scripts/buildkite"
fi

# ------------ Helpers ------------
require_file() {
  local f="${1:?missing path}"
  if [[ -r "$f" ]]; then
    return 0
  fi
  if [[ "$DRY_RUN" = "1" ]]; then
    echo "[DRY-RUN] would need file: $f (skipping strict check)" >&2
    return 0
  fi
  echo "Error: required file not found/readable: $f" >&2
  exit 1
}

file_exists() {
  [[ -e "$1" ]]
}

# Extract a quoted value from a CMake-style line: NAME "value"
get_field() {
  local fname="${1:?missing file}" field="${2:?missing field}"
  require_file "$fname"
  if ! file_exists "$fname"; then
    # DRY-RUN fallback
    echo ""
    return 0
  fi
  sed -n -E "s/.*${field}[[:space:]]*\"([^\"]+)\".*/\1/p" "$fname" | head -n1
}

get_version() {
  local vfile="${FEELPP_CMAKE_DIR}/feelpp.version.cmake"
  require_file "$vfile"

  local major minor micro prerelease
  major="$(get_field "$vfile" "VERSION_MAJOR")"
  minor="$(get_field "$vfile" "VERSION_MINOR")"
  micro="$(get_field "$vfile" "VERSION_MICRO")"
  prerelease="$(get_field "$vfile" "VERSION_PRERELEASE" || true)"

  if [[ -z "${major:-}" || -z "${minor:-}" || -z "${micro:-}" ]]; then
    if [[ "$DRY_RUN" = "1" ]]; then
      echo "[DRY-RUN] synthesizing FEELPP version as v0.0.0" >&2
      echo "v0.0.0"
      return 0
    fi
    echo "Error: could not parse version from $vfile" >&2
    exit 1
  fi

  if [[ -n "${prerelease:-}" ]]; then
    printf "v%s\n" "${major}.${minor}.${micro}${prerelease}"
  else
    printf "v%s\n" "${major}.${minor}.${micro}"
  fi
}

# Split "os:tag" -> echoes "<os> <tag>"
split_from() {
  local in="${1-}"
  if [[ -z "$in" ]]; then
    in="${DEFAULT_TARGET}"
  fi

  local os tag
  IFS=':' read -r os tag <<<"$in"

  if [[ -z "${os:-}" || -z "${tag:-}" ]]; then
    echo "Error: invalid FROM spec '${in}', expected 'os:version' (e.g. debian:13)" >&2
    exit 1
  fi

  printf '%s %s\n' "$os" "$tag"
}

# ------------ list.sh access (with DRY-RUN fallback) ------------
_list_sh() {
  echo "${FEELPP_SCRIPTS_DIR}/list.sh"
}

have_list_sh() {
  local list="$(_list_sh)"
  [[ -x "$list" ]]
}

# Returns first-column image name matching "<branch>-<version>-<os>-<tag>"
tag_from_target() {
  local from="${1-}" branch="${2-}" version="${3-}"

  if [[ -z "${branch:-}" || -z "${version:-}" ]]; then
    echo "Error: tag_from_target requires branch and version" >&2
    exit 1
  fi

  local fromos fromtag
  read -r fromos fromtag < <(split_from "$from")

  local list="$(_list_sh)"
  if have_list_sh; then
    "$list" "$branch" "$version" \
      | awk -v b="$branch" -v v="$version" -v o="$fromos" -v t="$fromtag" '
          BEGIN { pat="^" b "-" v "-" o "-" t "($|[[:space:]])" }
          $0 ~ pat { print $1; exit }
        '
    return 0
  fi

  if [[ "$DRY_RUN" = "1" ]]; then
    # Synthesize a plausible tag
    local synth="${branch}-${version}-${fromos}-${fromtag}"
    echo "[DRY-RUN] synthesizing tag_from_target -> ${synth}" >&2
    printf "%s\n" "$synth"
    return 0
  fi

  require_file "$list" # hard fail in non-dry-run
}

# Returns "<os>-<tag>"
tag_from_os() {
  local from="${1-}"
  local fromos fromtag
  read -r fromos fromtag < <(split_from "$from")
  printf "%s-%s\n" "$fromos" "$fromtag"
}

# Returns up to 3 extra tags (columns 2..4) for same match as tag_from_target
extratags_from_target() {
  local from="${1-}" branch="${2-}" version="${3-}"

  if [[ -z "${branch:-}" || -z "${version:-}" ]]; then
    echo "Error: extratags_from_target requires branch and version" >&2
    exit 1
  fi

  local fromos fromtag
  read -r fromos fromtag < <(split_from "$from")

  local list="$(_list_sh)"
  if have_list_sh; then
    "$list" "$branch" "$version" \
      | awk -v b="$branch" -v v="$version" -v o="$fromos" -v t="$fromtag" '
          BEGIN { pat="^" b "-" v "-" o "-" t "($|[[:space:]])" }
          $0 ~ pat {
            out="";
            for (i=2; i<=4 && i<=NF; i++) {
              if ($i != "") {
                if (out != "") out = out " ";
                out = out $i
              }
            }
            print out;
            exit
          }
        '
    return 0
  fi

  if [[ "$DRY_RUN" = "1" ]]; then
    echo "[DRY-RUN] synthesizing extratags_from_target -> (none)" >&2
    echo ""   # no extra tags synthesized by default
    return 0
  fi

  require_file "$list" # hard fail otherwise
}

# Combines a dockerfile template with a generated FROM line
# - Rewrites any explicit "FROM ghcr.io/feelpp/feelpp-env:*" inside the template
# - Rewrites any "ARG BASE=ghcr.io/feelpp/feelpp-env:*" too
# - If the template has no FROM at all, we prepend one
dockerfile_from() {
  local dockerfile="$1" from="$2"
  require_file "$dockerfile"

  # Does the template already reference our base image?
  if grep -Eq '^\s*FROM\s+ghcr\.io/feelpp/feelpp-env:' "$dockerfile" \
     || grep -Eq '^\s*ARG\s+BASE\s*=\s*ghcr\.io/feelpp/feelpp-env:' "$dockerfile"; then
    # Rewrite in-place stream (print to stdout)
    sed -E \
      -e "s|^(\s*FROM\s+)ghcr\.io/feelpp/feelpp-env:[^[:space:]]+|\1${from}|g" \
      -e "s|^(\s*ARG\s+BASE\s*=\s*)ghcr\.io/feelpp/feelpp-env:[^[:space:]]+|\1${from}|g" \
      "$dockerfile"
  else
    # If the template has no FROM lines, just prepend one
    if ! grep -Eq '^\s*FROM\s+' "$dockerfile"; then
      printf 'FROM %s\n' "$from"
      cat "$dockerfile"
    else
      # Template has FROMs but not feelpp-env ones; safest is to still prepend ours
      printf 'FROM %s\n' "$from"
      cat "$dockerfile"
    fi
  fi
}

# Export FEELPP_VERSION for callers
FEELPP_VERSION="$(get_version)"
export FEELPP_VERSION