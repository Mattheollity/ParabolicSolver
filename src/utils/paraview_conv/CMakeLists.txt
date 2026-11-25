add_executable(pvconv
    src/field_reader.cpp
    src/paraview_conv.cpp
    src/vtk_writer.cpp
)

target_include_directories(pvconv PRIVATE include)
target_compile_features(pvconv PRIVATE cxx_std_17)

set_target_properties(pvconv PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
