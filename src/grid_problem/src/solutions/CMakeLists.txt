cmake_minimum_required(VERSION 3.20)

add_library(
    solutions
    STATIC
    src/initial_conditions.cpp
    src/boundary_conditions.cpp
)

target_include_directories(solutions PUBLIC include)

target_link_libraries(solutions PUBLIC common)