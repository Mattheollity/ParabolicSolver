cmake_minimum_required(VERSION 3.20)

add_library(space_scheme STATIC
    src/space_scheme.cpp
    src/space_matrix_assembler.cpp
)

target_include_directories(space_scheme PUBLIC
    include
)

target_link_libraries(space_scheme PUBLIC
    common
    grid_problem
)