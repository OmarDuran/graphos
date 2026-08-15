#!/usr/bin/env bash
#
# Build the graphos portability stack -- camp, RAJA, Umpire, CHAI -- into one
# prefix, pinned to the coordinated LLNL 2024.07 release line so the four agree
# on camp and BLT. This is the Spack environment's specs (see spack.yaml)
# expressed as a plain source build, for machines that do not have Spack.
#
# ONE DEFINITION, TWO CONSUMERS. The CI workflow and the Dockerfile both call
# this script rather than each spelling the builds out, because a stack that is
# described twice is a stack that drifts: the container would keep passing on
# flags CI had already changed, and neither would be evidence about the other.
#
#   scripts/build_tpls.sh [PREFIX]
#
# PREFIX defaults to $TPL_PREFIX, then to ./tpl. Everything lands under it as a
# normal install tree, so a consumer needs only -DCMAKE_PREFIX_PATH=$PREFIX.
set -euo pipefail

PREFIX="${1:-${TPL_PREFIX:-$PWD/tpl}}"
VERSION="${TPL_VERSION:-2024.07.0}"
CAMP_VERSION="${CAMP_VERSION:-2024.07.0}"
JOBS="${JOBS:-$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu || echo 4)}"
SRC="${TPL_SRC:-$(mktemp -d)}"
# OpenMP on RAJA only: its for_policy becomes omp_parallel_for_exec, so every
# kernel-form op parallelizes with no code change. Umpire and CHAI keep it off
# -- their OpenMP variants gate internal threading the exec seams never use.
RAJA_OPENMP="${RAJA_OPENMP:-ON}"

mkdir -p "$PREFIX"
PREFIX="$(cd "$PREFIX" && pwd)"

echo "=== graphos TPL stack ${VERSION} -> ${PREFIX} (${JOBS} jobs) ==="

# A shallow clone WITH submodules: BLT is a submodule of each of these, and
# camp is a submodule of RAJA and Umpire. We build camp once ourselves and
# point the rest at it, which is what the Spack line does and what keeps a
# single camp in the prefix instead of three that merely happen to match.
clone() {  # clone <repo> <tag> <dir>
  [ -d "$3" ] || git clone --depth 1 --branch "$2" --recurse-submodules \
    --shallow-submodules "https://github.com/LLNL/$1.git" "$3"
}

common_flags=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
  -DCMAKE_PREFIX_PATH="$PREFIX"
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  -DENABLE_TESTS=OFF
  -DENABLE_EXAMPLES=OFF
  -DENABLE_BENCHMARKS=OFF
  -DENABLE_DOCS=OFF
)

build() {  # build <src-dir> [extra cmake args...]
  local src="$1"; shift
  cmake -S "$src" -B "$src/build" "${common_flags[@]}" "$@"
  cmake --build "$src/build" -j "$JOBS"
  cmake --install "$src/build"
}

# ---- camp: the shared type layer the other three agree on -------------------
clone camp "v${CAMP_VERSION}" "$SRC/camp"
build "$SRC/camp"

# ---- RAJA: kernel execution ------------------------------------------------
clone RAJA "v${VERSION}" "$SRC/raja"
build "$SRC/raja" \
  -DENABLE_OPENMP="$RAJA_OPENMP" \
  -DRAJA_ENABLE_OPENMP="$RAJA_OPENMP" \
  -DRAJA_ENABLE_TESTS=OFF \
  -DRAJA_ENABLE_EXAMPLES=OFF \
  -DRAJA_ENABLE_BENCHMARKS=OFF \
  -Dcamp_DIR="$PREFIX/lib/cmake/camp"

# ---- Umpire: memory pools --------------------------------------------------
# ~fmt_header_only, as the Spack line requires: header-only fmt exports
# FMT_HEADER_ONLY to consumers and inlines its implementation into every
# translation unit, which breaks under C++20 consteval with recent compilers on
# the fmt version this Umpire pins. Linking the compiled fmt avoids the path.
clone Umpire "v${VERSION}" "$SRC/umpire"
build "$SRC/umpire" \
  -DENABLE_OPENMP=OFF \
  -DUMPIRE_ENABLE_TOOLS=OFF \
  -DUMPIRE_ENABLE_FMT_HEADER_ONLY=OFF \
  -Dcamp_DIR="$PREFIX/lib/cmake/camp"

# ---- CHAI: managed arrays --------------------------------------------------
clone CHAI "v${VERSION}" "$SRC/chai"
build "$SRC/chai" \
  -DENABLE_OPENMP=OFF \
  -DCHAI_ENABLE_TESTS=OFF \
  -DCHAI_ENABLE_EXAMPLES=OFF \
  -DRAJA_DIR="$PREFIX/lib/cmake/raja" \
  -Dumpire_DIR="$PREFIX/lib/cmake/umpire" \
  -Dcamp_DIR="$PREFIX/lib/cmake/camp"

echo "=== installed into ${PREFIX} ==="
ls "$PREFIX/lib/cmake"
