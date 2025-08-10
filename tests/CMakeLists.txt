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

add_executable(test_pml test_pml.cpp)
target_link_libraries(test_pml PRIVATE vectorem)
add_test(NAME test_pml COMMAND test_pml)
set_tests_properties(test_pml PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} LABELS pml)

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
