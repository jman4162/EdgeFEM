cmake_minimum_required(VERSION 3.20)
project(VectorEM VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 3.4 REQUIRED)

option(VECTOREM_BUILD_SCALAR "Build scalar prototype" ON)
option(VECTOREM_BUILD_VECTOR "Build vector Maxwell path" ON)
option(VECTOREM_PYTHON "Build Python bindings" OFF)

# When building via scikit-build-core, SKBUILD is defined
if(SKBUILD)
  set(VECTOREM_PYTHON ON CACHE BOOL "" FORCE)
endif()

add_subdirectory(src)
add_subdirectory(tests)
if (VECTOREM_PYTHON)
  add_subdirectory(python)
endif()

if (VECTOREM_BUILD_SCALAR)
  add_test(NAME smoke_scalar_demo COMMAND vectorem_scalar_demo ${PROJECT_SOURCE_DIR}/examples/cube_cavity.msh 1)
  set_tests_properties(smoke_scalar_demo PROPERTIES LABELS smoke)
endif()
