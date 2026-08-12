# Installing graphos

graphos is a **header-only C++20 library**. There is nothing to compile to
*use* it — installation copies headers and a CMake package so other projects
can `find_package(graphos)`. The optional hardware-portability stack (RAJA,
Umpire, CHAI) is only needed if you enable it.

Requirements: **CMake ≥ 3.20** and a **C++20 compiler** (Apple Clang, Clang,
or GCC). Nothing else for the default build.

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

Pick any prefix you like. Without `--prefix` the install goes to the
configure-time `CMAKE_INSTALL_PREFIX` (default `/usr/local`, which needs
`sudo`).

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

graphos is header-only — add the include directory and C++20:

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
| `GRAPHOS_FETCH_TPL` | `OFF` | Download+build RAJA/Umpire via FetchContent instead of `find_package` (development only; skips the install export, and CHAI is not fetchable this way) |
| `CMAKE_INSTALL_PREFIX` | `/usr/local` | Install destination (or use `cmake --install build --prefix …`) |

Each `GRAPHOS_ENABLE_*` option locates its library with `find_package`, so
the TPLs must be installed and discoverable (see the Spack section). The
installed `graphosConfig.cmake` remembers which options were on and
re-finds those dependencies for consumers automatically.

---

## 4. The portability stack via Spack

The repository ships a Spack environment ([spack.yaml](spack.yaml)) pinning
RAJA, Umpire, and CHAI to the coordinated **2024.07** LLNL release line so
the three agree on their shared dependencies (camp, BLT).

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

This builds RAJA, Umpire, and CHAI (plus camp/BLT) and links them into a
single **view** at `.spack-env/view/` inside the repo (git-ignored). First
build takes a while; afterwards the environment is cached.

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

Consumers of a portability-enabled install must also be able to find the
TPLs — the simplest way is to configure them with the same
`CMAKE_PREFIX_PATH` (view path first, then the graphos prefix):

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

---

## 5. Troubleshooting

**`find_package(graphos)` fails** — pass the install prefix:
`-DCMAKE_PREFIX_PATH=$HOME/opt/graphos` when configuring the consumer.

**`find_package(RAJA/umpire/chai)` fails with a `GRAPHOS_ENABLE_*` option** —
the TPLs are not discoverable. Activate the Spack environment or pass
`-DCMAKE_PREFIX_PATH=<repo>/.spack-env/view`.

**Two Spack installations fighting each other** (e.g. one from Homebrew and
a `~/spack` clone): they share `~/.spack` config but speak different config
dialects, which surfaces as parse errors like
`a single spec was requested, but parsed more than one`. Pick one — check
`which spack` and `spack --version`, and invoke the one you mean by full
path (`~/spack/bin/spack …`) or fix your `PATH`.

**Bootstrap error `No module named 'clingo.ast'`** — the binary-cached
concretizer does not match your Python. The most reliable fix is installing
clingo directly into the Python Spack runs on (bypasses bootstrap entirely):

```bash
$(spack python -c "import sys; print(sys.executable)") -m pip install clingo
```

**macOS: builds fail instantly with an empty log** (`make: no such file`) —
Apple's `/usr/bin/make` is the ancient 3.81 and Homebrew installs GNU make
as `gmake` only, while Spack invokes `<prefix>/bin/make`. Install
`brew install make`, create a shim providing the `make` name, and register
it:

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

**Install export was skipped** — you configured with `GRAPHOS_FETCH_TPL=ON`.
That mode is for development; install builds must find the TPLs via
`find_package` (the Spack view).
