include(FetchContent)
# PINNED. `stable` is a moving tag: it advanced to pybind11 3.x, whose
# stricter overload deduction broke this file with no change on our side.
# A fetched dependency that can move is a build that can break without a commit.
FetchContent_Declare(pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.13.6)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_core bindings.cpp)
target_link_libraries(_core PRIVATE graphos::graphos)
set_target_properties(_core PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/graphos)

if(BUILD_TESTING)
  add_test(NAME python.test_graphos
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_graphos.py)
  set_tests_properties(python.test_graphos PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME python.test_networkx
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_networkx.py)
  set_tests_properties(python.test_networkx PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}")
endif()
