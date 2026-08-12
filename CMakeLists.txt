cmake_minimum_required(VERSION 3.20)
project(graphos VERSION 0.1.0 LANGUAGES CXX)

add_library(graphos INTERFACE)
add_library(graphos::graphos ALIAS graphos)
target_include_directories(graphos INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(graphos INTERFACE cxx_std_20)

# Hardware portability seams. Off by default so the host fallback builds with
# nothing but a C++20 compiler; enabling swaps in the DOE stack.
option(GRAPHOS_ENABLE_RAJA "Dispatch exec::forall/scan through RAJA policies" OFF)
option(GRAPHOS_ENABLE_UMPIRE "Allocate exec::Buffer from Umpire memory pools" OFF)
option(GRAPHOS_ENABLE_CHAI "Store frozen complexes in CHAI managed arrays" OFF)
option(GRAPHOS_FETCH_TPL "FetchContent third-party libs instead of find_package" OFF)

if(GRAPHOS_ENABLE_RAJA)
  if(GRAPHOS_FETCH_TPL)
    include(FetchContent)
    set(RAJA_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(RAJA_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(RAJA_ENABLE_EXERCISES OFF CACHE BOOL "" FORCE)
    set(RAJA_ENABLE_BENCHMARKS OFF CACHE BOOL "" FORCE)
    set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(raja
      GIT_REPOSITORY https://github.com/LLNL/RAJA.git
      GIT_TAG v2024.07.0)
    FetchContent_MakeAvailable(raja)
  else()
    find_package(RAJA REQUIRED)
  endif()
  target_link_libraries(graphos INTERFACE RAJA)
  target_compile_definitions(graphos INTERFACE GRAPHOS_HAVE_RAJA)
endif()

if(GRAPHOS_ENABLE_UMPIRE)
  if(GRAPHOS_FETCH_TPL)
    include(FetchContent)
    FetchContent_Declare(umpire
      GIT_REPOSITORY https://github.com/LLNL/Umpire.git
      GIT_TAG v2024.07.0)
    FetchContent_MakeAvailable(umpire)
  else()
    find_package(umpire REQUIRED)
  endif()
  target_link_libraries(graphos INTERFACE umpire)
  target_compile_definitions(graphos INTERFACE GRAPHOS_HAVE_UMPIRE)
endif()

if(GRAPHOS_ENABLE_CHAI)
  # CHAI coordinates versions with RAJA/Umpire; use installed packages
  # (e.g. via spack), not GRAPHOS_FETCH_TPL
  find_package(chai REQUIRED)
  target_link_libraries(graphos INTERFACE chai)
  target_compile_definitions(graphos INTERFACE GRAPHOS_HAVE_CHAI)
endif()

include(CTest)
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  add_subdirectory(tests)
endif()

# ---- installation ----------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(GRAPHOS_FETCH_TPL AND (GRAPHOS_ENABLE_RAJA OR GRAPHOS_ENABLE_UMPIRE))
  # FetchContent'd TPLs have no install export to reference; fetch builds are
  # a development convenience, install against find_package'd TPLs instead
  message(STATUS "graphos: GRAPHOS_FETCH_TPL build - skipping CMake package export")
else()
  install(TARGETS graphos EXPORT graphosTargets)
  install(EXPORT graphosTargets
          FILE graphosTargets.cmake
          NAMESPACE graphos::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/graphos)
  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/graphosConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/graphosConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/graphos)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/graphosConfigVersion.cmake
    COMPATIBILITY SameMinorVersion
    ARCH_INDEPENDENT)
  install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/graphosConfig.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/graphosConfigVersion.cmake
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/graphos)
endif()
