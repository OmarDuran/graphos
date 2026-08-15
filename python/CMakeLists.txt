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
  # pybind11 exports the interpreter as Python_EXECUTABLE under FindPython and
  # as PYTHON_EXECUTABLE in classic mode (the pinned 2.13.6 uses the latter).
  # Taking only one of them leaves the test command without an interpreter, so
  # ctest invokes the .py file itself and reports BAD_COMMAND.
  set(GRAPHOS_PYTHON "${Python_EXECUTABLE}")
  if(NOT GRAPHOS_PYTHON)
    set(GRAPHOS_PYTHON "${PYTHON_EXECUTABLE}")
  endif()
  if(NOT GRAPHOS_PYTHON)
    message(FATAL_ERROR "graphos: python tests requested but no interpreter was resolved")
  endif()

  add_test(NAME python.test_graphos
           COMMAND ${GRAPHOS_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_graphos.py)
  set_tests_properties(python.test_graphos PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME python.test_networkx
           COMMAND ${GRAPHOS_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_networkx.py)
  set_tests_properties(python.test_networkx PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}")
endif()
