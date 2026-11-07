cmake_minimum_required(VERSION 3.10)
project(PoissonDecomposition CXX)



set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

add_subdirectory(libs/eigen)
add_subdirectory(libs/Catch2)
add_subdirectory(libs/eigen-expokit)

add_subdirectory(src/grid)
add_subdirectory(src/space_scheme)
add_subdirectory(src/time_scheme)
add_subdirectory(src/utils)
add_subdirectory(src/paraview_conv)

add_executable(compt
    src/main.cpp
)

target_link_libraries(
    compt
    grid
    poisson
    time_scheme
    utils
)

