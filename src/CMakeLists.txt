add_executable(vectorem_scalar_demo scalar_demo.cpp)

target_include_directories(vectorem_scalar_demo PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(vectorem_scalar_demo PRIVATE Eigen3::Eigen)
