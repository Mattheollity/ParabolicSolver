cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver_file_io CXX)



set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

add_library(file_io
    src/basic_io/basic_io.cpp
    src/basic_io/eigen_io.cpp
    src/modules_io/configs_io.cpp
    src/modules_io/geometry_io.cpp
    src/modules_io/solutions_io.cpp
    src/modules_io/space_matrix_assembler_io.cpp
    src/modules_io/space_scheme_io.cpp
)

target_include_directories(file_io PUBLIC include/basic_io)
target_include_directories(file_io PUBLIC include/modules_io)
target_link_libraries(file_io PUBLIC time_scheme space_scheme grid)
target_link_libraries(file_io PUBLIC common)
target_link_libraries(file_io PRIVATE stdc++fs)
target_link_libraries(file_io PRIVATE solutions)