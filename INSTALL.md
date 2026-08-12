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
concretizer does not match your Python. Bootstrap from source instead:

```bash
spack bootstrap disable github-actions-v2
```

```bash
spack bootstrap now
```

**Install export was skipped** — you configured with `GRAPHOS_FETCH_TPL=ON`.
That mode is for development; install builds must find the TPLs via
`find_package` (the Spack view).
