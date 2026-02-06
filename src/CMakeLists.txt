add_library(edgefem
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

set_target_properties(edgefem PROPERTIES POSITION_INDEPENDENT_CODE ON)

if (EDGEFEM_BUILD_SCALAR)
  target_sources(edgefem PRIVATE fem_scalar.cpp)
endif()

if (EDGEFEM_BUILD_VECTOR)
  target_sources(edgefem PRIVATE edge_basis.cpp assemble_maxwell.cpp)
  add_executable(edgefem_waveguide_demo waveguide_demo.cpp)
  target_link_libraries(edgefem_waveguide_demo PRIVATE edgefem)
  add_executable(edgefem_sweep_demo sweep_demo.cpp)
  target_link_libraries(edgefem_sweep_demo PRIVATE edgefem)
endif()

target_include_directories(edgefem PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(edgefem PUBLIC Eigen3::Eigen)

if (EDGEFEM_BUILD_SCALAR)
  add_executable(edgefem_scalar_demo main_scalar_demo.cpp)
  target_link_libraries(edgefem_scalar_demo PRIVATE edgefem)
endif()
