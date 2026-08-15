# unit tests mirror include/graphos/: one test file per header, one ctest
# entry per file, named <folder>.<header>
set(GRAPHOS_UNIT_TESTS
  core/test_types
  core/test_build
  core/test_complex
  core/test_coboundary
  core/test_frozen
  core/test_incidence
  core/test_marker
  core/test_union_find
  exec/test_forall
  exec/test_memory
  exec/test_array
  ops/test_disjoint_union
  ops/test_quotient
  ops/test_pushout
  ops/test_star_deletion
  ops/test_cut
  ops/test_subcomplex
  queries/test_facets
  queries/test_components
  ops/test_product
  ops/test_lift_identifications
  ops/test_subdivision
  ops/test_orient
  queries/test_homology
  queries/test_manifold
  ops/test_agglomerate
  ops/test_join
  queries/test_neighborhood
  queries/test_adjacency
  queries/test_amalgamation
  ops/test_collapse
  ops/test_coarsen
  ops/test_replace
  property/test_properties
  property/test_scaling
  property/test_performance
  property/test_multilevel
  property/test_compliance
)

foreach(t IN LISTS GRAPHOS_UNIT_TESTS)
  string(REPLACE "/" "." test_name ${t})
  add_executable(${test_name} ${t}.cpp)
  target_link_libraries(${test_name} PRIVATE graphos::graphos)
  target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ${CMAKE_SOURCE_DIR}/bench)
  if(NOT MSVC)
    target_compile_options(${test_name} PRIVATE -Wall -Wextra)
  endif()
  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
