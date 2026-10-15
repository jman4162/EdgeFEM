# FindUMFPACK
# -----------
# Finds the SuiteSparse UMFPACK library and header.
#
# SuiteSparse 7 ships UMFPACKConfig.cmake, but its SuiteSparse_config
# dependency requires OpenMP for C, which AppleClang does not provide
# without extra flags. This module looks for the header and the library
# directly instead (Homebrew: /opt/homebrew/include/suitesparse/umfpack.h;
# Debian/Ubuntu libsuitesparse-dev: /usr/include/suitesparse/umfpack.h).
# The shared libumfpack records its own dependencies (AMD, CHOLMOD,
# SuiteSparse_config, BLAS), so linking it alone resolves the umfpack_*
# calls Eigen's UmfPackSupport makes.
#
# Result variables: UMFPACK_FOUND, UMFPACK_INCLUDE_DIR, UMFPACK_LIBRARY,
# UMFPACK_VERSION. Imported target: UMFPACK::UMFPACK.
#
# Set CMAKE_DISABLE_FIND_PACKAGE_UMFPACK=ON to configure as if UMFPACK
# were not installed.

find_path(UMFPACK_INCLUDE_DIR umfpack.h PATH_SUFFIXES suitesparse)
find_library(UMFPACK_LIBRARY NAMES umfpack)

if (UMFPACK_INCLUDE_DIR AND EXISTS "${UMFPACK_INCLUDE_DIR}/umfpack.h")
  file(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" _umf_version_lines
       REGEX "#define UMFPACK_(MAIN|SUB|SUBSUB)_VERSION")
  foreach (_part MAIN SUB SUBSUB)
    string(REGEX REPLACE ".*#define UMFPACK_${_part}_VERSION +([0-9]+).*"
           "\\1" _umf_${_part} "${_umf_version_lines}")
  endforeach()
  set(UMFPACK_VERSION "${_umf_MAIN}.${_umf_SUB}.${_umf_SUBSUB}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK
  REQUIRED_VARS UMFPACK_LIBRARY UMFPACK_INCLUDE_DIR
  VERSION_VAR UMFPACK_VERSION)

if (UMFPACK_FOUND AND NOT TARGET UMFPACK::UMFPACK)
  add_library(UMFPACK::UMFPACK UNKNOWN IMPORTED)
  set_target_properties(UMFPACK::UMFPACK PROPERTIES
    IMPORTED_LOCATION "${UMFPACK_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIR}")
endif()

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY)
