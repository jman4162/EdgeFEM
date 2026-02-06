find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(pyedgefem pyedgefem.cpp)
target_link_libraries(pyedgefem PRIVATE edgefem)
target_include_directories(pyedgefem PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Install target for scikit-build-core wheel builds
install(TARGETS pyedgefem LIBRARY DESTINATION .)
