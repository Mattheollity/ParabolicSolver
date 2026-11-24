cmake_minimum_required(VERSION 3.20)

add_library(
    interactions
    STATIC
    src/grid.cpp
    src/boundaries_interaction.cpp
    src/self_subdomain_interaction.cpp
)

target_include_directories(interactions PUBLIC include)

target_link_libraries(interactions PUBLIC common)
target_link_libraries(interactions PUBLIC geometry)
target_link_libraries(interactions PUBLIC solutions)