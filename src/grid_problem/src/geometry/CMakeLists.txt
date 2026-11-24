cmake_minimum_required(VERSION 3.20)

add_library(
    geometry
    STATIC
    src/base_geometry/vertex.cpp
    src/base_geometry/dim_iterator.cpp
    src/geometry_creation/domain_geometry_creation.cpp
    src/geometry_creation/subdomain_creation.cpp
    src/geometry.cpp
)

target_include_directories(geometry PUBLIC include)

target_link_libraries(geometry PUBLIC common)