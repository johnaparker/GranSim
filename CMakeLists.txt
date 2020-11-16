cmake_minimum_required(VERSION 3.1)
project(cpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-Werror -Og")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
add_subdirectory(pybind11)

set(SOURCES ${SOURCES} 
    src/gran2d.cpp
    src/gran3d.cpp
    src/geometry_2d.cpp
    src/geometry_3d.cpp
)

set(BINDINGS
    src/main_pyb.cpp
    src/gran2d_pyb.cpp
    src/gran3d_pyb.cpp
)

add_library(cpp MODULE "${SOURCES}" "${BINDINGS}")
target_link_libraries(cpp PRIVATE pybind11::module OpenMP::OpenMP_CXX Eigen3::Eigen)

set_target_properties(cpp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")
