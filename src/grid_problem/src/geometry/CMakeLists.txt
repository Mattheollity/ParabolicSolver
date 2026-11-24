cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver_geometry CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    geometry
    STATIC
    src/base_geometry/vertex.cpp
    src/base_geometry/dim_iterator.cpp
    src/geometry_creation/domain_geometry_creation.cpp
    src/geometry_creation/subdomain_creation.cpp
    src/geometry.cpp
)

target_include_directories(
    geometry
    PUBLIC
    include
)

target_link_libraries(geometry PUBLIC common)