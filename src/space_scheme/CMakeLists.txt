cmake_minimum_required(VERSION 3.10)
project(PoissonApprox CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(
    space_scheme
    STATIC
    src/matrix_assembler.cpp
    src/space_scheme.cpp
)
target_link_libraries(space_scheme PUBLIC grid)

target_include_directories(
    space_scheme
    PUBLIC
    include
)



# add_executable(
#     poisson_test
#     tests/poisson.cpp
# )
# target_link_libraries(poisson_test PUBLIC poisson)
# target_link_libraries(poisson_test PRIVATE Catch2::Catch2WithMain)