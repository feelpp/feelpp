#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <core|electric|mor|heat|solid|hdg>" >&2
  exit 2
fi

lane="$1"
jobs="${JOBS:-20}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"

case "$lane" in
  core)
    configure_preset="python-core-ci"
    build_preset="python-core-ci"
    test_preset="python-core-ci"
    build_dir="build/python-core-ci"
    install_prefix="install/python-core-ci"
    build_targets=(pyfeelpp)
    smoke_modules=(feelpp.core)
    ;;
  electric)
    configure_preset="python-electric-ci"
    build_preset="python-electric-ci"
    test_preset="python-electric-ci"
    build_dir="build/python-electric-ci"
    install_prefix="install/python-electric-ci"
    build_targets=(pyfeelpp pyfeelpptoolboxes)
    smoke_modules=(feelpp.core feelpp.toolboxes.core feelpp.toolboxes.electric)
    ;;
  mor)
    configure_preset="python-mor-ci"
    build_preset="python-mor-ci"
    test_preset="python-mor-ci"
    build_dir="build/python-mor-ci"
    install_prefix="install/python-mor-ci"
    build_targets=(pyfeelpp pyfeelpptoolboxes pyfeelppmor)
    smoke_modules=(feelpp.core feelpp.toolboxes.core feelpp.toolboxes.electric feelpp.mor feelpp.mor.online feelpp.mor.reducedbasis.reducedbasis)
    ;;
  heat)
    configure_preset="python-heat-ci"
    build_preset="python-heat-ci"
    test_preset="python-heat-ci"
    build_dir="build/python-heat-ci"
    install_prefix="install/python-heat-ci"
    build_targets=(pyfeelpp pyfeelpptoolboxes)
    smoke_modules=(feelpp.core feelpp.toolboxes.core feelpp.toolboxes.heat feelpp.core.interpolation)
    ;;
  solid)
    configure_preset="python-solid-ci"
    build_preset="python-solid-ci"
    test_preset="python-solid-ci"
    build_dir="build/python-solid-ci"
    install_prefix="install/python-solid-ci"
    build_targets=(pyfeelpp pyfeelpptoolboxes)
    smoke_modules=(feelpp.core feelpp.toolboxes.core feelpp.toolboxes.solid)
    ;;
  hdg)
    configure_preset="python-hdg-ci"
    build_preset="python-hdg-ci"
    test_preset="python-hdg-ci"
    build_dir="build/python-hdg-ci"
    install_prefix="install/python-hdg-ci"
    build_targets=(pyfeelpp pyfeelpptoolboxes)
    smoke_modules=(feelpp.core feelpp.toolboxes.core feelpp.toolboxes.hdg)
    ;;
  *)
    echo "unknown lane: $lane" >&2
    exit 2
    ;;
esac

rm -rf "$build_dir" "$install_prefix"

cd "$repo_root"
cmake --preset "$configure_preset"
cmake --build --preset "$build_preset" -j "$jobs" --target "${build_targets[@]}"
ctest --preset "$test_preset"
cmake --install "$build_dir" --component Libs
cmake --install "$build_dir" --component fmt_core
cmake --install "$build_dir" --component Python

build_dir_abs="$repo_root/$build_dir"
install_prefix_abs="$repo_root/$install_prefix"

python_module_path="$(sed -n 's/^FEELPP_PYTHON_MODULE_PATH:STRING=//p' "$build_dir_abs/CMakeCache.txt" | head -n 1)"
if [[ -z "$python_module_path" ]]; then
  echo "failed to read FEELPP_PYTHON_MODULE_PATH from $build_dir_abs/CMakeCache.txt" >&2
  exit 1
fi

python_executable="$(sed -n 's/^Python3_EXECUTABLE:FILEPATH=//p' "$build_dir_abs/CMakeCache.txt" | head -n 1)"
if [[ -z "$python_executable" ]]; then
  echo "failed to read Python3_EXECUTABLE from $build_dir_abs/CMakeCache.txt" >&2
  exit 1
fi

python_site_packages="$install_prefix_abs/$python_module_path"
smoke_cmd=("$python_executable" "$repo_root/scripts/ci/python_install_smoke.py" --prefix "$install_prefix_abs")
for module_name in "${smoke_modules[@]}"; do
  smoke_cmd+=(--module "$module_name")
done

(
  cd /tmp
  env \
    HOME=/tmp \
    PYTHONNOUSERSITE=1 \
    PYTHONPATH="$python_site_packages" \
    LD_LIBRARY_PATH="$install_prefix_abs/lib" \
    "${smoke_cmd[@]}"
)
