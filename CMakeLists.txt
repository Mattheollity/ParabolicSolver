cmake_minimum_required(VERSION 3.20)
project(ParabolicSolver CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

add_subdirectory(libs/eigen)

add_library(common INTERFACE)

target_include_directories(common INTERFACE include/common)

target_link_libraries(common INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(common INTERFACE Eigen3::Eigen)

add_subdirectory(src/grid_problem/src/solutions)
add_subdirectory(src/grid_problem/src/geometry)