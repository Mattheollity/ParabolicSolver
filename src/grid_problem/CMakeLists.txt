cmake_minimum_required(VERSION 3.20)

add_subdirectory(src/geometry)
add_subdirectory(src/solutions)

add_library(grid_problem STATIC
    src/grid_problem.cpp
)

target_include_directories(grid_problem PUBLIC
    include
)

target_link_libraries(grid_problem PUBLIC
    common
    geometry
    solutions
)