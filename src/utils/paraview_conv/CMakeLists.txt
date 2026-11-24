cmake_minimum_required(VERSION 3.10)
project(ParaviewConverter CXX)

add_executable(
    pvconv
    src/field_reader.cpp
    src/paraview_conv.cpp
    src/vtk_writer.cpp
)

target_include_directories(pvconv PUBLIC include)
target_link_libraries(pvconv PUBLIC file_io)
target_link_libraries(pvconv PUBLIC common)

set_target_properties(pvconv PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)