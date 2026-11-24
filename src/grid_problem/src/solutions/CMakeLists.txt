cmake_minimum_required(VERSION 3.20)

add_library(
    solutions
    STATIC
    src/solutions.cpp
)

target_include_directories(solutions PUBLIC include)

target_link_libraries(solutions PUBLIC common)