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

add_subdirectory(src/utils/mpi_utils)
add_subdirectory(src/grid_problem)
add_subdirectory(src/computation/space_scheme)
add_subdirectory(src/computation/time_scheme)
add_subdirectory(src/utils/file_io)
add_subdirectory(src/utils/paraview_conv)

add_library(parabolic_solver STATIC
    src/parabolic_solver.cpp
)

target_include_directories(parabolic_solver PUBLIC include)

target_link_libraries(parabolic_solver PUBLIC
    common
    file_io
)



add_executable(ParabolicSolver
    src/main.cpp
)

target_link_libraries(ParabolicSolver PRIVATE
    parabolic_solver
)