cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver CXX)



set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

add_library(common INTERFACE)
target_include_directories(common INTERFACE include/common)
target_link_libraries(common INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(common INTERFACE Eigen3::Eigen)

add_subdirectory(libs/eigen)
add_subdirectory(libs/Catch2)
add_subdirectory(libs/eigen-expokit)

add_subdirectory(src/solutions)
add_subdirectory(src/grid)
add_subdirectory(src/space_scheme)
add_subdirectory(src/mpi_utils)
add_subdirectory(src/time_scheme)
add_subdirectory(src/file_io)
add_subdirectory(src/paraview_conv)

add_library(
    parabolic_solver
    src/parabolic_solver.cpp
)

target_include_directories(parabolic_solver PRIVATE include)

target_link_libraries(
    parabolic_solver
    PRIVATE
    common
    file_io
    mpi_utils
    grid
    space_scheme
    time_scheme
    paraview_conv
)

add_executable(ParabolicSolver
    src/main.cpp
)

target_link_libraries(ParabolicSolver PRIVATE parabolic_solver)