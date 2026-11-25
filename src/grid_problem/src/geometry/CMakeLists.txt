cmake_minimum_required(VERSION 3.20)

add_library(
    geometry
    STATIC
    base_geometry/src/vertex.cpp
    base_geometry/src/geometry.cpp
    base_geometry/src/interactions.cpp
    cartesian_uniform/src/cartesian_uniform.cpp
    cartesian_uniform/src/dim_iterator.cpp
    cartesian_uniform/src/cartesian_uniform_creation.cpp
    cartesian_uniform/src/cartesian_uniform_subdomain_creation.cpp
)

target_include_directories(
    geometry
    PUBLIC
    base_geometry/include
    cartesian_uniform/include
)

target_link_libraries(geometry PUBLIC common)