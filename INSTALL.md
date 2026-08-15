# Installing graphos

graphos is a **header-only C++20 library**. Nothing has to be compiled to use
it: installing copies the headers and a CMake package so other projects can
`find_package(graphos)`. The portability stack (RAJA, Umpire, CHAI) is needed
only if you enable it.

Requirements: **CMake ≥ 3.20** and a **C++20 compiler** (Apple Clang, Clang or
GCC). Nothing else for the default build.

---

## 1. Quick install (no dependencies)

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

```bash
cmake --build build -j 8
```

```bash
ctest --test-dir build
```

```bash
cmake --install build --prefix ~/opt/graphos
```

Any prefix works. Without `--prefix` the install goes to the configure-time
`CMAKE_INSTALL_PREFIX`, default `/usr/local`, which needs `sudo`.

### What gets installed

```
<prefix>/include/graphos/            all headers (graphos.hpp is the umbrella)
<prefix>/lib/cmake/graphos/          graphosConfig.cmake, version + targets files
```

---

## 2. Using the installed package

### From CMake (recommended)

```cmake
find_package(graphos CONFIG REQUIRED)
target_link_libraries(myapp PRIVATE graphos::graphos)
```

If the prefix is not a system location, tell CMake where it is when
configuring *your* project:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=$HOME/opt/graphos
```

### Without CMake

Header-only, so the include directory and C++20 are enough:

```bash
c++ -std=c++20 -I$HOME/opt/graphos/include my_program.cpp -o my_program
```

---

## 3. Configure options

| Option | Default | Effect |
|---|---|---|
| `GRAPHOS_ENABLE_RAJA` | `OFF` | `exec::forall`/scans dispatch through RAJA policies |
| `GRAPHOS_ENABLE_UMPIRE` | `OFF` | `exec::Buffer` allocates from Umpire memory pools |
| `GRAPHOS_ENABLE_CHAI` | `OFF` | Frozen-complex storage (`exec::Array`) uses `chai::ManagedArray` |
| `GRAPHOS_FETCH_TPL` | `OFF` | Download and build RAJA/Umpire with FetchContent instead of `find_package`. Development only: it skips the install export, and CHAI cannot be fetched this way |
| `GRAPHOS_BUILD_BENCH` | `ON` | Build `graphos_bench` |
| `GRAPHOS_BUILD_PYTHON` | `OFF` | Build the `graphos._core` bindings into `python/graphos/` (fetches pybind11, pinned) |
| `GRAPHOS_SANITIZE` | *(empty)* | Sanitizers for tests and bench, e.g. `address,undefined`. Top-level builds only |
| `BUILD_TESTING` | `ON` | Standard CTest switch; `OFF` skips the test tree |
| `CMAKE_INSTALL_PREFIX` | `/usr/local` | Install destination, or use `cmake --install build --prefix …` |

Each `GRAPHOS_ENABLE_*` option locates its library with `find_package`, so the
TPLs must be installed and discoverable — see §4 and §5. The installed
`graphosConfig.cmake` records which options were on and re-finds those
dependencies for consumers.

---

## 4. The portability stack via Spack

The repository ships a Spack environment ([spack.yaml](spack.yaml)) pinning
RAJA, Umpire and CHAI to the coordinated **2024.07** LLNL release line, so the
three agree on their shared dependencies (camp, BLT). If you would rather not
use Spack, §5 builds the same line from source with one script.

### One-time Spack setup

If you do not have Spack yet:

```bash
git clone --depth=2 https://github.com/spack/spack.git ~/spack
```

```bash
source ~/spack/share/spack/setup-env.sh
```

(Put the `source` line in your shell profile to make `spack` permanent.)

### Build the TPLs

From the graphos repository root:

```bash
spack env activate .
```

```bash
spack install
```

This builds RAJA, Umpire and CHAI (with camp and BLT) and links them into a
single **view** at `.spack-env/view/`, git-ignored. The first build is slow;
afterwards the environment is cached.

### Build graphos against the stack

```bash
cmake -S . -B build-portable -DCMAKE_BUILD_TYPE=Release \
  -DGRAPHOS_ENABLE_RAJA=ON -DGRAPHOS_ENABLE_UMPIRE=ON -DGRAPHOS_ENABLE_CHAI=ON \
  -DCMAKE_PREFIX_PATH=$PWD/.spack-env/view
```

```bash
cmake --build build-portable -j 8
```

```bash
ctest --test-dir build-portable
```

```bash
cmake --install build-portable --prefix ~/opt/graphos
```

A consumer of a portability-enabled install must also find the TPLs. The
simplest way is the same `CMAKE_PREFIX_PATH`, view first, then the graphos
prefix:

```bash
cmake -S . -B build "-DCMAKE_PREFIX_PATH=/path/to/graphos-repo/.spack-env/view;$HOME/opt/graphos"
```

### GPU variants

Edit the three specs in `spack.yaml` to carry the matching variant, e.g. for
CUDA on an H100:

```yaml
  - raja@2024.07.0 +cuda cuda_arch=90
  - umpire@2024.07.0 +cuda cuda_arch=90
  - chai@2024.07.0 +cuda cuda_arch=90
```

(`+rocm amdgpu_target=gfx90a` for AMD.) Then `spack install` again.

Note that no device execution policy exists yet — building the stack with
`+cuda` prepares the dependencies, but `exec::forall` still dispatches to host
policies. See the Roadmap in [README.md](README.md).

---

## 5. The portability stack from source (no Spack)

[scripts/build_tpls.sh](scripts/build_tpls.sh) builds the same 2024.07 line —
camp, RAJA, Umpire, CHAI — into one prefix, using only git and CMake:

```bash
scripts/build_tpls.sh /opt/graphos-tpl
```

camp is built first and the other three are pointed at it, so the prefix holds
one camp rather than three that happen to match. `JOBS`, `TPL_VERSION` and
`RAJA_OPENMP` override the defaults.

Then configure graphos against it exactly as with the Spack view:

```bash
cmake -S . -B build-portable -DCMAKE_BUILD_TYPE=Release \
  -DGRAPHOS_ENABLE_RAJA=ON -DGRAPHOS_ENABLE_UMPIRE=ON -DGRAPHOS_ENABLE_CHAI=ON \
  -DCMAKE_PREFIX_PATH=/opt/graphos-tpl
```

The CI workflow and the container call this same script, so all three build
the stack one way. Editing the script invalidates the CI cache, which is why
the cache key is its hash.

---

## 6. The container

[Dockerfile](Dockerfile) runs the whole thing from a bare image, in the four
stages the CI workflow uses:

```bash
docker build -t graphos .                    # every stage
docker build -t graphos --target tested .    # stop after ctest
docker build -t graphos --build-arg ENABLE_TPL=OFF .   # the no-TPL path
```

The stages are separate so a failure names its phase — `tpls` is a
third-party problem, `library` is ours, `tested` is a behaviour change, and
`consumer` means the install exports are wrong although everything built and
passed. The final image carries the install at `/opt/graphos` with
`CMAKE_PREFIX_PATH` already set.

---

## 7. Troubleshooting

**`find_package(graphos)` fails** — pass the install prefix:
`-DCMAKE_PREFIX_PATH=$HOME/opt/graphos` when configuring the consumer.

**`find_package(RAJA/umpire/chai)` fails with a `GRAPHOS_ENABLE_*` option** —
the TPLs are not discoverable. Activate the Spack environment or pass
`-DCMAKE_PREFIX_PATH=<repo>/.spack-env/view`.

**Two Spack installations in conflict** (say one from Homebrew and a
`~/spack` clone): they share `~/.spack` but read different config dialects,
which surfaces as `a single spec was requested, but parsed more than one`.
Choose one — check `which spack` and `spack --version`, then invoke it by full
path or fix `PATH`.

**Bootstrap error `No module named 'clingo.ast'`** — the binary-cached
concretizer does not match your Python. The most reliable fix is installing
clingo directly into the Python Spack runs on (bypasses bootstrap entirely):

```bash
$(spack python -c "import sys; print(sys.executable)") -m pip install clingo
```

**macOS: builds fail immediately with an empty log** (`make: no such file`) —
Apple ships make 3.81 and Homebrew installs GNU make as `gmake` only, while
Spack invokes `<prefix>/bin/make`. Install `brew install make`, then shim the
`make` name and register it:

```bash
mkdir -p ~/.spack/shims/bin && ln -sf /opt/homebrew/bin/gmake ~/.spack/shims/bin/make
```

then in `~/.spack/packages.yaml`:

```yaml
packages:
  gmake:
    externals:
    - spec: gmake@4.4.1
      prefix: ~/.spack/shims
    buildable: false
```

**macOS: TPL builds fail with `unsupported option '-mcpu='` or
`-arch x86_64` appearing on an arm64 machine** — Spack is running on an
Intel (Rosetta) Python, typically from an x86_64 miniconda
(`file $(spack python -c "import sys; print(sys.executable)")` to check).
Every universal binary the build spawns (shell, clang, CMake) then prefers
its x86_64 slice, so CMake targets x86_64 while Spack emits arm64 tuning
flags. Point Spack at a native arm64 Python (and give it clingo):

```bash
/opt/homebrew/bin/python3 -m pip install --break-system-packages clingo
```

```bash
export SPACK_PYTHON=/opt/homebrew/bin/python3
```

(Add the export to your shell profile so it sticks.)

**macOS: `consteval function ... is not a constant expression` errors from
`fmt/format-inl.h` when building with Umpire/CHAI** — Umpire's default
`+fmt_header_only` exports `FMT_HEADER_ONLY` to consumers, inlining fmt's
implementation into every C++20 translation unit, which breaks with recent
AppleClang and the fmt ≤ 11.0 this Umpire line pins. The bundled
[spack.yaml](spack.yaml) sets `umpire ~fmt_header_only` for this reason;
keep that variant if you edit the environment.

**Install export was skipped** — the build was configured with
`GRAPHOS_FETCH_TPL=ON`, which is a development mode. An install build must
find the TPLs through `find_package`: the Spack view (§4) or the source prefix
(§5).

**The Python bindings fail to build with overload errors** — pybind11 is
pinned in `python/CMakeLists.txt` for this reason. `GIT_TAG stable` is a
moving tag and advanced to 3.x, whose stricter overload deduction breaks the
bindings with no change on this side. Keep the pin.
