add_library(vectorem
  mesh_gmsh.cpp
  fem_scalar.cpp
  solver.cpp
  bc.cpp
)

target_include_directories(vectorem PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(vectorem PUBLIC Eigen3::Eigen)

add_executable(vectorem_scalar_demo main_scalar_demo.cpp)

target_link_libraries(vectorem_scalar_demo PRIVATE vectorem)

