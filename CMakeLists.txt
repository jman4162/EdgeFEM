cmake_minimum_required(VERSION 3.20)

# The version lives in python/edgefem/__init__.py (__version__); pyproject.toml
# and this file both read it from there.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/python/edgefem/__init__.py"
     _edgefem_version_line REGEX "^__version__ = \"[0-9]+\\.[0-9]+\\.[0-9]+\"")
string(REGEX REPLACE "^__version__ = \"([0-9.]+)\".*" "\\1"
       EDGEFEM_VERSION "${_edgefem_version_line}")
if (NOT EDGEFEM_VERSION)
  message(FATAL_ERROR "no __version__ = \"X.Y.Z\" line in python/edgefem/__init__.py")
endif()

project(EdgeFEM VERSION ${EDGEFEM_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 3.4 REQUIRED)

option(EDGEFEM_BUILD_SCALAR "Build scalar prototype" ON)
option(EDGEFEM_BUILD_VECTOR "Build vector Maxwell path" ON)
option(EDGEFEM_PYTHON "Build Python bindings" OFF)
option(EDGEFEM_BUILD_DEBUG_TOOLS "Build one-off diagnostic executables in tests/" OFF)
option(EDGEFEM_USE_SUITESPARSE
       "Use SuiteSparse UMFPACK (GPL-2.0-or-later) for direct solves" OFF)

# Direct solves use Eigen SparseLU unless EDGEFEM_USE_SUITESPARSE=ON.
# UMFPACK is GPL-2.0-or-later, so it is linked only on request, never
# because it happens to be installed; a request without UMFPACK is an
# error rather than a silent fallback.
set(EDGEFEM_HAVE_UMFPACK OFF)
if (EDGEFEM_USE_SUITESPARSE)
  find_package(UMFPACK MODULE)
  if (NOT UMFPACK_FOUND)
    message(FATAL_ERROR
      "EDGEFEM_USE_SUITESPARSE=ON but SuiteSparse UMFPACK was not found "
      "(looked for umfpack.h, also under include/suitesparse, and "
      "libumfpack). Install it (macOS: brew install suite-sparse; "
      "Debian/Ubuntu: apt install libsuitesparse-dev), point "
      "CMAKE_PREFIX_PATH at it, or configure with "
      "-DEDGEFEM_USE_SUITESPARSE=OFF to use Eigen SparseLU.")
  endif()
  set(EDGEFEM_HAVE_UMFPACK ON)
  message(STATUS "EdgeFEM: direct solves use SuiteSparse UMFPACK "
                 "${UMFPACK_VERSION} (GPL-2.0-or-later)")
else()
  message(STATUS "EdgeFEM: direct solves use Eigen SparseLU "
                 "(EDGEFEM_USE_SUITESPARSE=OFF)")
endif()

# When building via scikit-build-core, SKBUILD is defined
if(SKBUILD)
  set(EDGEFEM_PYTHON ON CACHE BOOL "" FORCE)
endif()

add_subdirectory(src)
add_subdirectory(tests)
if (EDGEFEM_PYTHON)
  add_subdirectory(python)
endif()

if (EDGEFEM_BUILD_SCALAR)
  add_test(NAME smoke_scalar_demo COMMAND edgefem_scalar_demo ${PROJECT_SOURCE_DIR}/examples/cube_cavity.msh 1)
  set_tests_properties(smoke_scalar_demo PROPERTIES LABELS smoke)
endif()
