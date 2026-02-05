add_executable(test_edge_indexing test_edge_indexing.cpp)
target_link_libraries(test_edge_indexing PRIVATE vectorem)
add_test(NAME test_edge_indexing COMMAND test_edge_indexing)
set_tests_properties(test_edge_indexing PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_maxwell test_maxwell.cpp)
target_link_libraries(test_maxwell PRIVATE vectorem)
add_test(NAME test_maxwell COMMAND test_maxwell)
set_tests_properties(test_maxwell PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_ports test_ports.cpp)
target_link_libraries(test_ports PRIVATE vectorem)
add_test(NAME test_ports COMMAND test_ports)
set_tests_properties(test_ports PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_waveguide_modal test_waveguide_modal.cpp)
target_link_libraries(test_waveguide_modal PRIVATE vectorem)
add_test(NAME test_waveguide_modal COMMAND test_waveguide_modal)
set_tests_properties(test_waveguide_modal PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_pml test_pml.cpp)
target_link_libraries(test_pml PRIVATE vectorem)
add_test(NAME test_pml COMMAND test_pml)
set_tests_properties(test_pml PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS pml)

add_executable(test_abc test_abc.cpp)
target_link_libraries(test_abc PRIVATE vectorem)
add_test(NAME test_abc COMMAND test_abc)
set_tests_properties(test_abc PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS abc)

add_executable(test_sweep test_sweep.cpp)
target_link_libraries(test_sweep PRIVATE vectorem)
add_test(NAME test_sweep COMMAND test_sweep)
set_tests_properties(test_sweep PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS smoke)

add_executable(test_vector_fit test_vector_fit.cpp)
target_link_libraries(test_vector_fit PRIVATE vectorem)
add_test(NAME test_vector_fit COMMAND test_vector_fit)
set_tests_properties(test_vector_fit PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS smoke)

add_executable(test_ntf test_ntf.cpp)
target_link_libraries(test_ntf PRIVATE vectorem)
add_test(NAME test_ntf COMMAND test_ntf)
set_tests_properties(test_ntf PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_solver_convergence test_solver_convergence.cpp)
target_link_libraries(test_solver_convergence PRIVATE vectorem)
add_test(NAME test_solver_convergence COMMAND test_solver_convergence)
set_tests_properties(test_solver_convergence PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS smoke)

add_executable(test_mesh_quality test_mesh_quality.cpp)
target_link_libraries(test_mesh_quality PRIVATE vectorem)
add_test(NAME test_mesh_quality COMMAND test_mesh_quality)
set_tests_properties(test_mesh_quality PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_touchstone_nport test_touchstone_nport.cpp)
target_link_libraries(test_touchstone_nport PRIVATE vectorem)
add_test(NAME test_touchstone_nport COMMAND test_touchstone_nport)
set_tests_properties(test_touchstone_nport PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_coupling test_coupling.cpp)
target_link_libraries(test_coupling PRIVATE vectorem)
add_test(NAME test_coupling COMMAND test_coupling)
set_tests_properties(test_coupling PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_periodic test_periodic.cpp)
target_link_libraries(test_periodic PRIVATE vectorem)
add_test(NAME test_periodic COMMAND test_periodic)
set_tests_properties(test_periodic PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(debug_waveguide debug_waveguide.cpp)
target_link_libraries(debug_waveguide PRIVATE vectorem)

add_executable(benchmark_2d_waveguide benchmark_2d_waveguide.cpp)
target_link_libraries(benchmark_2d_waveguide PRIVATE vectorem)
add_test(NAME benchmark_2d_waveguide COMMAND benchmark_2d_waveguide)
set_tests_properties(benchmark_2d_waveguide PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(debug_fem_core debug_fem_core.cpp)
target_link_libraries(debug_fem_core PRIVATE vectorem)

add_executable(debug_port_weights debug_port_weights.cpp)
target_link_libraries(debug_port_weights PRIVATE vectorem)

add_executable(test_frequency_sweep test_frequency_sweep.cpp)
target_link_libraries(test_frequency_sweep PRIVATE vectorem)

add_executable(test_normalized_ports test_normalized_ports.cpp)
target_link_libraries(test_normalized_ports PRIVATE vectorem)

add_executable(debug_port_normals debug_port_normals.cpp)
target_link_libraries(debug_port_normals PRIVATE vectorem)

add_executable(test_flipped_port test_flipped_port.cpp)
target_link_libraries(test_flipped_port PRIVATE vectorem)

add_executable(debug_matrix_cond debug_matrix_cond.cpp)
target_link_libraries(debug_matrix_cond PRIVATE vectorem)

add_executable(test_weight_normalization test_weight_normalization.cpp)
target_link_libraries(test_weight_normalization PRIVATE vectorem)

add_executable(test_waveguide_modes test_waveguide_modes.cpp)
target_link_libraries(test_waveguide_modes PRIVATE vectorem)

add_executable(debug_port_coupling debug_port_coupling.cpp)
target_link_libraries(debug_port_coupling PRIVATE vectorem)

add_executable(debug_mode_alignment debug_mode_alignment.cpp)
target_link_libraries(debug_mode_alignment PRIVATE vectorem)

add_executable(debug_3d_mode_structure debug_3d_mode_structure.cpp)
target_link_libraries(debug_3d_mode_structure PRIVATE vectorem)

add_executable(debug_solution_structure debug_solution_structure.cpp)
target_link_libraries(debug_solution_structure PRIVATE vectorem)

add_executable(test_2d_port_eigensolve test_2d_port_eigensolve.cpp)
target_link_libraries(test_2d_port_eigensolve PRIVATE vectorem)

add_executable(test_waveport_with_abc test_waveport_with_abc.cpp)
target_link_libraries(test_waveport_with_abc PRIVATE vectorem)

add_executable(debug_3d_vs_2d_pattern debug_3d_vs_2d_pattern.cpp)
target_link_libraries(debug_3d_vs_2d_pattern PRIVATE vectorem)

add_executable(test_waveguide_with_pml test_waveguide_with_pml.cpp)
target_link_libraries(test_waveguide_with_pml PRIVATE vectorem)

add_executable(test_3d_eigenvector_excitation test_3d_eigenvector_excitation.cpp)
target_link_libraries(test_3d_eigenvector_excitation PRIVATE vectorem)

add_executable(test_eigenvector_waveport test_eigenvector_waveport.cpp)
target_link_libraries(test_eigenvector_waveport PRIVATE vectorem)

add_executable(test_sparam_accuracy test_sparam_accuracy.cpp)
target_link_libraries(test_sparam_accuracy PRIVATE vectorem)
add_test(NAME test_sparam_accuracy COMMAND test_sparam_accuracy)
set_tests_properties(test_sparam_accuracy PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_direct_mode_excitation test_direct_mode_excitation.cpp)
target_link_libraries(test_direct_mode_excitation PRIVATE vectorem)

add_executable(test_direct_mode_tuning test_direct_mode_tuning.cpp)
target_link_libraries(test_direct_mode_tuning PRIVATE vectorem)

add_executable(test_eigenmode_sparams test_eigenmode_sparams.cpp)
target_link_libraries(test_eigenmode_sparams PRIVATE vectorem)
add_test(NAME test_eigenmode_sparams COMMAND test_eigenmode_sparams)
set_tests_properties(test_eigenmode_sparams PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS smoke)

add_executable(benchmark_wr90 benchmark_wr90.cpp)
target_link_libraries(benchmark_wr90 PRIVATE vectorem)
add_test(NAME benchmark_wr90 COMMAND benchmark_wr90)
set_tests_properties(benchmark_wr90 PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS benchmark)

if (VECTOREM_PYTHON)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME python_smoke
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python_smoke.py)
  set_tests_properties(python_smoke PROPERTIES
                       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                       ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python:${PROJECT_SOURCE_DIR}"
                       LABELS smoke)

  add_test(NAME python_sdk_demo
           COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/examples/sdk_demo.py)
  set_tests_properties(python_sdk_demo PROPERTIES
                       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                       ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python:${PROJECT_SOURCE_DIR}"
                       LABELS smoke)
endif()
