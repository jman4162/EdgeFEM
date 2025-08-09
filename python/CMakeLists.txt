find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(pyvectorem pyvectorem.cpp)
target_link_libraries(pyvectorem PRIVATE vectorem)
target_include_directories(pyvectorem PRIVATE ${PROJECT_SOURCE_DIR}/include)
