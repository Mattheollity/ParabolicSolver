cmake_minimum_required(VERSION 3.20)

add_subdirectory(geometry)
add_subdirectory(interactions)
add_subdirectory(solutions)

add_library(
    grid_problem
    STATIC
    src/grid_problem.cpp
)

target_include_directories(grid_problem PUBLIC include)

target_link_libraries(grid_problem PUBLIC common)
target_link_libraries(grid_problem PUBLIC geometry)
target_link_libraries(grid_problem PUBLIC interactions)
target_link_libraries(grid_problem PUBLIC solutions)