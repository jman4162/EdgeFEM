add_library(vectorem
  mesh_gmsh.cpp
  mesh_quality.cpp
  solver.cpp
  bc.cpp
  coupling.cpp
  periodic.cpp
  ports/port_eigensolve.cpp
  ports/wave_port.cpp
  io/touchstone.cpp
  io/json_export.cpp
  io/array_export.cpp
  io/vtk_export.cpp
  array/embedded_pattern.cpp
  array/active_impedance.cpp
  array/coupling_extract.cpp
  post/ntf.cpp
  sweep.cpp
  mor/vector_fit.cpp
  materials/dispersive.cpp
)

set_target_properties(vectorem PROPERTIES POSITION_INDEPENDENT_CODE ON)

if (VECTOREM_BUILD_SCALAR)
  target_sources(vectorem PRIVATE fem_scalar.cpp)
endif()

if (VECTOREM_BUILD_VECTOR)
  target_sources(vectorem PRIVATE edge_basis.cpp assemble_maxwell.cpp)
  add_executable(vectorem_waveguide_demo waveguide_demo.cpp)
  target_link_libraries(vectorem_waveguide_demo PRIVATE vectorem)
  add_executable(vectorem_sweep_demo sweep_demo.cpp)
  target_link_libraries(vectorem_sweep_demo PRIVATE vectorem)
endif()

target_include_directories(vectorem PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(vectorem PUBLIC Eigen3::Eigen)

if (VECTOREM_BUILD_SCALAR)
  add_executable(vectorem_scalar_demo main_scalar_demo.cpp)
  target_link_libraries(vectorem_scalar_demo PRIVATE vectorem)
endif()

