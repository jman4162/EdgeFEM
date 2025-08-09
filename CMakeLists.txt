cmake_minimum_required(VERSION 3.20)
project(VectorEM VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 3.4 REQUIRED)

add_subdirectory(src)

add_test(NAME smoke_scalar_demo COMMAND vectorem_scalar_demo ${PROJECT_SOURCE_DIR}/examples/cube_cavity.msh 1)
set_tests_properties(smoke_scalar_demo PROPERTIES LABELS smoke)
