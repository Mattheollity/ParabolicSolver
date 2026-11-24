cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver_grid CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src/geometry)

add_library(
    grid
    STATIC
    src/grid.cpp
)

target_include_directories(
    grid
    PUBLIC
    include
)

target_link_libraries(grid PUBLIC geometry)
target_link_libraries(grid PUBLIC common)
target_link_libraries(grid PUBLIC solutions)