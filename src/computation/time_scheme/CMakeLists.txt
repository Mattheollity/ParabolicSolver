add_library(time_scheme STATIC
    src/time_scheme.cpp
    src/time_step_subdomain.cpp
    src/time_step_boundary.cpp
    src/time_scheme_assembler.cpp
    src/time_scheme_assembler_mpi.cpp
)

target_include_directories(time_scheme PUBLIC include)

target_link_libraries(time_scheme PUBLIC
    common
    space_scheme
    mpi_utils
)