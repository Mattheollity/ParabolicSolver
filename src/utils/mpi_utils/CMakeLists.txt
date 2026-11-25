find_package(MPI REQUIRED COMPONENTS CXX)

set(MPI_CXX_COMPILER "/usr/bin/mpicxx.mpich"
    CACHE FILEPATH "MPICH C++ wrapper")
set(MPI_HOME "/usr/lib/x86_64-linux-gnu/mpich")

add_library(mpi_utils STATIC
    src/mpi_utils.cpp
)

target_include_directories(mpi_utils PUBLIC include)

target_link_libraries(mpi_utils PUBLIC
    common
    MPI::MPI_CXX
)