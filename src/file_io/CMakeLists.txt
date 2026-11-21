cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver_file_io CXX)



set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

add_library(file_io
    src/basic_io.cpp
    src/configs_io.cpp
    src/eigen_io.cpp
    src/geometry_io.cpp
)

target_include_directories(file_io PUBLIC include)
target_link_libraries(file_io PUBLIC time_scheme space_scheme grid)
target_link_libraries(file_io PUBLIC common)
target_link_libraries(file_io PRIVATE stdc++fs)
target_link_libraries(file_io PRIVATE solutions)