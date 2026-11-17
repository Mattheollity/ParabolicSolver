cmake_minimum_required(VERSION 3.10)
project(NonOverlapDecomposeGrid CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    grid
    STATIC
    src/grid.cpp
    src/grid_interactions.cpp
)
target_link_libraries(grid PUBLIC geometry)
target_include_directories(
    grid
    PUBLIC
    include
)



# add_executable(
#     grid_test
#     tests/grid.cpp
# )

# target_link_libraries(grid_test PUBLIC grid)
# target_link_libraries(grid_test PRIVATE Catch2::Catch2WithMain)

# set_target_properties(grid_test PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
# )