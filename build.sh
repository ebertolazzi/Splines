#!/usr/bin/env bash

# This script is used to build the project using CMake.
# Usage:
#   ./build.sh [option] [build_type] [-p install_prefix]
# Where:
#   option: "configure", "build", "test", "install", "package" (default: build)
#   build_type: "Debug", "Release", "RelWithDebInfo", "MinSizeRel" (default: Release)
#   install_prefix: path to install the project (default: /usr/local)
#
# Build directory is "build" by default. You can change it by setting the BUILD_DIR environment variable.
# Generator can be specified by setting the CMAKE_GENERATOR environment variable (default: Ninja).

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./build.sh [option] [build_type] [-p install_prefix]

Where:
  option:         configure, build, test, install, package (default: build)
  build_type:     Debug, Release, RelWithDebInfo, MinSizeRel (default: Release)
  install_prefix: path to install the project (default: /usr/local)

Environment:
  BUILD_DIR        Build directory (default: build)
  CMAKE_GENERATOR  CMake generator (default: Ninja)

Examples:
  ./build.sh configure Debug
  ./build.sh build Release
  ./build.sh test
  ./build.sh install Release -p "$HOME/.local"
  ./build.sh package
EOF
}

die() {
  echo "build.sh: $*" >&2
  echo >&2
  usage >&2
  exit 1
}

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_root"

option="build"
build_type="Release"
install_prefix="/usr/local"
option_seen=0
build_type_seen=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -p|--prefix)
      [[ $# -ge 2 ]] || die "missing argument for $1"
      install_prefix="$2"
      shift 2
      ;;
    configure|build|test|install|package)
      [[ $option_seen -eq 0 ]] || die "option specified more than once"
      option="$1"
      option_seen=1
      shift
      ;;
    Debug|Release|RelWithDebInfo|MinSizeRel)
      [[ $build_type_seen -eq 0 ]] || die "build_type specified more than once"
      build_type="$1"
      build_type_seen=1
      shift
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

build_dir="${BUILD_DIR:-build}"
generator="${CMAKE_GENERATOR:-Ninja}"

configure() {
  cmake \
    -S . \
    -B "$build_dir" \
    -G "$generator" \
    -DCMAKE_BUILD_TYPE="$build_type" \
    -DCMAKE_INSTALL_PREFIX="$install_prefix"
}

ensure_configured() {
  if [[ ! -f "$build_dir/CMakeCache.txt" ]]; then
    configure
  fi
}

case "$option" in
  configure)
    configure
    ;;
  build)
    ensure_configured
    cmake --build "$build_dir" --config "$build_type"
    ;;
  test)
    ensure_configured
    cmake --build "$build_dir" --config "$build_type"
    ctest --test-dir "$build_dir" --build-config "$build_type"
    ;;
  install)
    ensure_configured
    cmake --build "$build_dir" --config "$build_type"
    cmake --install "$build_dir" --config "$build_type" --prefix "$install_prefix"
    ;;
  package)
    ensure_configured
    cmake --build "$build_dir" --config "$build_type" --target package
    ;;
esac
