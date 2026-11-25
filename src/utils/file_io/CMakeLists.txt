add_library(file_io STATIC
    src/basic_io/basic_io.cpp
    src/basic_io/eigen_io.cpp

    src/modules_io/configs_io.cpp
    src/modules_io/solutions_io.cpp
    src/modules_io/geometry_io.cpp
    src/modules_io/space_scheme_io.cpp
    src/modules_io/space_matrix_assembler_io.cpp
    src/modules_io/time_scheme_io.cpp
)

target_include_directories(file_io PUBLIC
    include/basic_io
    include/modules_io
)

target_link_libraries(file_io PUBLIC
    common
    time_scheme
)