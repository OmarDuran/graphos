# graphos in a container, built in the four stages the CI workflow runs, in the
# same order and with the same scripts. The point of the image is not to ship a
# binary -- graphos is header-only -- but to be an INSTALLATION that can be
# checked: a machine with nothing on it, a stack built from source, a package
# installed, and a consumer that finds it with find_package.
#
#   docker build -t graphos .                   # the full stack (slow, cached)
#   docker build -t graphos --target tested .   # stop after the test stage
#   docker build -t graphos --build-arg ENABLE_TPL=OFF .   # the no-TPL path
#
# The stages are separate so a failure names its phase. A break in `tpls` is a
# third-party problem; a break in `library` is ours; a break in `tested` is a
# behaviour change; a break in `consumer` means the install exports are wrong
# even though everything built and passed.

ARG UBUNTU_VERSION=24.04

# ---- phase 1: the third-party stack -----------------------------------------
FROM ubuntu:${UBUNTU_VERSION} AS tpls

ARG ENABLE_TPL=ON
ENV DEBIAN_FRONTEND=noninteractive TPL_PREFIX=/opt/tpl

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# only the script, so editing a header does not invalidate the stack layer --
# this is the layer worth caching, and it is the expensive one
COPY scripts/build_tpls.sh /usr/local/bin/build_tpls.sh
RUN chmod +x /usr/local/bin/build_tpls.sh \
    && if [ "$ENABLE_TPL" = "ON" ]; then \
         TPL_SRC=/tmp/tpl-src /usr/local/bin/build_tpls.sh "$TPL_PREFIX" \
         && rm -rf /tmp/tpl-src; \
       else \
         mkdir -p "$TPL_PREFIX" && echo "TPL stack skipped (ENABLE_TPL=OFF)"; \
       fi

# ---- phase 2: the library ---------------------------------------------------
FROM tpls AS library

ARG ENABLE_TPL=ON
ARG BUILD_TYPE=Release
WORKDIR /src
COPY . .

# CHAI is reached through find_package only, so it rides on the prefix rather
# than on GRAPHOS_FETCH_TPL -- which is why the stack is built first.
RUN cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
      -DCMAKE_INSTALL_PREFIX=/opt/graphos \
      -DCMAKE_PREFIX_PATH=${TPL_PREFIX} \
      -DGRAPHOS_ENABLE_RAJA=${ENABLE_TPL} \
      -DGRAPHOS_ENABLE_UMPIRE=${ENABLE_TPL} \
      -DGRAPHOS_ENABLE_CHAI=${ENABLE_TPL} \
    && cmake --build build -j "$(nproc)" \
    && cmake --install build

# ---- phase 3: the tests -----------------------------------------------------
FROM library AS tested
RUN ctest --test-dir build --output-on-failure \
    && ./build/bench/graphos_bench 8

# ---- phase 4: the install is usable by a consumer ---------------------------
#
# Everything above can pass while the installed package is unusable: a missing
# export, an absolute path baked into a config file, a dependency the target
# does not carry. So the last stage throws the build tree away and compiles a
# fresh consumer against the INSTALL, which is the only thing a user ever sees.
FROM tested AS consumer
RUN rm -rf /src/build \
    && mkdir -p /tmp/consumer && cd /tmp/consumer \
    && printf '%s\n' \
       'cmake_minimum_required(VERSION 3.20)' \
       'project(consumer LANGUAGES CXX)' \
       'find_package(graphos CONFIG REQUIRED)' \
       'add_executable(consumer main.cpp)' \
       'target_link_libraries(consumer PRIVATE graphos::graphos)' > CMakeLists.txt \
    && printf '%s\n' \
       '#include <cstdio>' \
       '#include "graphos/graphos.hpp"' \
       'int main() {' \
       '  graphos::Complex c(2);' \
       '  c.attach_vertices(3);' \
       '  c.attach_cell(1, {0, 1}, {-1, 1});' \
       '  std::printf("consumer: dim=%d edges=%lld\\n", c.dim(), (long long)c.count(1));' \
       '  return c.count(1) == 1 ? 0 : 1;' \
       '}' > main.cpp \
    && cmake -S . -B build -DCMAKE_PREFIX_PATH="/opt/graphos;${TPL_PREFIX}" \
    && cmake --build build \
    && ./build/consumer

# ---- the shipped image ------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION} AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      libstdc++6 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=consumer /opt/graphos /opt/graphos
COPY --from=consumer /opt/tpl /opt/tpl
ENV CMAKE_PREFIX_PATH=/opt/graphos:/opt/tpl
LABEL org.opencontainers.image.title="graphos" \
      org.opencontainers.image.description="A metric-free computational topology engine"
