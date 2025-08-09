add_executable(test_edge_indexing test_edge_indexing.cpp)
target_link_libraries(test_edge_indexing PRIVATE vectorem)
add_test(NAME test_edge_indexing COMMAND test_edge_indexing)
set_tests_properties(test_edge_indexing PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(test_maxwell test_maxwell.cpp)
target_link_libraries(test_maxwell PRIVATE vectorem)
add_test(NAME test_maxwell COMMAND test_maxwell)
set_tests_properties(test_maxwell PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
