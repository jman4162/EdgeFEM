# Minimal Eigen3 finder
# Tries Config mode first, then falls back to header search
find_package(Eigen3 3.4 QUIET CONFIG)
if(Eigen3_FOUND)
    return()
endif()

find_path(EIGEN3_INCLUDE_DIR Eigen/Dense PATH_SUFFIXES eigen3)
if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 3.4 or later not found")
endif()

add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
