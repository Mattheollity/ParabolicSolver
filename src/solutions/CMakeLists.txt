cmake_minimum_required(VERSION 3.10)
project(ParabolicSolver_solutions CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

add_library(solutions ./solutions.cpp)
target_include_directories(solutions PUBLIC include)
target_link_libraries(solutions PUBLIC common)