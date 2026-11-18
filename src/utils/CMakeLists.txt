cmake_minimum_required(VERSION 3.10)
project(NonOverlapDecomposeUtils CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

set(MPI_CXX_COMPILER "/usr/bin/mpicxx.mpich"
    CACHE FILEPATH "MPICH C++ wrapper")
set(MPI_HOME "/usr/lib/x86_64-linux-gnu/mpich"
    CACHE PATH "MPICH installation prefix")

find_package(MPI REQUIRED COMPONENTS CXX)

add_library(
    utils
    STATIC
    src/file_io.cpp
    src/problem_conditions.cpp
    src/mpi_utils.cpp
)

target_include_directories(
    utils
    PUBLIC
    include
)

target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(utils PUBLIC Eigen3::Eigen)
target_link_libraries(utils PUBLIC eigen-expokit::eexpokit)
target_link_libraries(utils PUBLIC stdc++fs)
target_link_libraries(utils PUBLIC MPI::MPI_CXX)
