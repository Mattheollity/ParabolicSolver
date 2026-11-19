cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver_geometry CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    geometry
    STATIC
    src/vertex.cpp
    src/dim_iterator.cpp
    src/geometry_io.cpp
    src/geometry.cpp
)

target_include_directories(
    geometry
    PUBLIC
    include
)

target_link_libraries(mpi_utils INTERFACE common)