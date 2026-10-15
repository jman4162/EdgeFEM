# Minimal Eigen3 finder
# Tries Config mode first, then falls back to header search.
#
# Eigen 3.4.x and 5.x both work. Each Eigen3ConfigVersion.cmake accepts only
# its own major version (a 3.4.0 config rejects a "3.4...5" range because
# the range's upper end is major 5), so each major version is asked for in
# turn.
find_package(Eigen3 3.4 QUIET CONFIG)
if(NOT Eigen3_FOUND)
    find_package(Eigen3 5 QUIET CONFIG)
endif()
if(Eigen3_FOUND)
    return()
endif()

find_path(EIGEN3_INCLUDE_DIR Eigen/Dense PATH_SUFFIXES eigen3)
if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 3.4 or later not found")
endif()

add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
