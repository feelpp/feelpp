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

# --- Arch detection -----------------------------------------------------------
# Returns buildx platforms for a given TARGET like "ubuntu:24.04"
arches_for_target() {
  local target="$1"
  case "$target" in
    ubuntu:24.04|debian:13)
      printf "linux/amd64,linux/arm64"
      ;;
    *)
      printf "linux/amd64"
      ;;
  esac
}

# Optional: make a nice OCI description using detected arches
description_for() {
  local image="$1" target="$2"
  local arches; arches="$(arches_for_target "$target")"
  printf "Feel++ container for %s (%s) — platforms: %s" "$image" "$target" "$arches"
}

# dockerfile_from <template> <from_base> [description]
dockerfile_from() {
  local dockerfile="$1" from="$2" desc="${3:-}"
  require_file "$dockerfile"

  local label_line=""
  if [[ -n "$desc" ]]; then
    label_line="LABEL org.opencontainers.image.description=\"${desc}\""
  fi

  _emit_from_and_label() {
    printf 'FROM %s\n' "$from"
    [[ -n "$label_line" ]] && printf '%s\n' "$label_line"
  }

  if ! grep -Eq '^\s*FROM\s+' "$dockerfile"; then
    _emit_from_and_label
    cat "$dockerfile"
    return
  fi

  # Replace first FROM if it refers to feelpp-env; otherwise prepend ours.
  if grep -Eq '^\s*FROM\s+ghcr\.io/feelpp/feelpp-env:' "$dockerfile"; then
    awk -v from="$from" -v label="$label_line" '
      BEGIN{done=0}
      /^[[:space:]]*FROM[[:space:]]+ghcr\.io\/feelpp\/feelpp-env:/ && !done {
        print "FROM " from
        if (label != "") print label
        done=1; next
      }
      { print }
    ' "$dockerfile"
  else
    _emit_from_and_label
    cat "$dockerfile"
  fi
}

dockerfile_from_with_label() {
  local template="$1" base="$2" description="$3"
  dockerfile_from "$template" "$base" \
    | sed "1a LABEL org.opencontainers.image.description=\"$description\""
}

# Export FEELPP_VERSION for callers
FEELPP_VERSION="$(get_version)"
export FEELPP_VERSION