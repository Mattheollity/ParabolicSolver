#
# CMake for eigen-expokit
#
cmake_minimum_required (VERSION 3.21)

# Default build mode
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
endif ()

project (eigen-expokit VERSION 0.1.0 LANGUAGES CXX)

option (EXPOKIT_BUILD_FORTRAN "Build the reference Fortran expokit and samples" OFF)
option (EXPOKIT_BUILD_TESTING "Build the example programs and the tests" ${PROJECT_IS_TOP_LEVEL})
option (EXPOKIT_INSTALL "Generate the install and package-config rules" ${PROJECT_IS_TOP_LEVEL})

# googletest requires C++17, the library itself only needs C++14
if (NOT CMAKE_CXX_STANDARD)
  set (CMAKE_CXX_STANDARD 17)
endif ()
set (CMAKE_CXX_STANDARD_REQUIRED ON)

include (GNUInstallDirs)

#
# The library: headers only
#
add_library (eexpokit INTERFACE)
add_library (eigen-expokit::eexpokit ALIAS eexpokit)

target_include_directories (eexpokit INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features (eexpokit INTERFACE cxx_std_14)

if (EXPOKIT_INSTALL)
  include (CMakePackageConfigHelpers)

  set (EXPOKIT_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/eigen-expokit")

  install (TARGETS eexpokit EXPORT eigen-expokit-targets)
  install (DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/eexpokit/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/eexpokit"
    FILES_MATCHING PATTERN "*.hpp"
    PATTERN "example.hpp" EXCLUDE)
  install (EXPORT eigen-expokit-targets
    FILE eigen-expokit-targets.cmake
    NAMESPACE eigen-expokit::
    DESTINATION "${EXPOKIT_CMAKEDIR}")

  configure_package_config_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/eigen-expokit-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/eigen-expokit-config.cmake"
    INSTALL_DESTINATION "${EXPOKIT_CMAKEDIR}")
  # 0.x minor releases are not promised to be compatible with each other
  write_basic_package_version_file (
    "${CMAKE_CURRENT_BINARY_DIR}/eigen-expokit-config-version.cmake"
    COMPATIBILITY SameMinorVersion
    ARCH_INDEPENDENT)
  install (FILES
    "${CMAKE_CURRENT_BINARY_DIR}/eigen-expokit-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/eigen-expokit-config-version.cmake"
    DESTINATION "${EXPOKIT_CMAKEDIR}")
endif ()

#set (SANITIZER ${SANITIZER} -O2 -fsanitize=address -fno-omit-frame-pointer)
#set (SANLIB -lasan)
#set (SANITIZER ${SANITIZER} -fsanitize=leak)
#target_link_libraries (example -lleak)
#set (SANITIZER "${SANITIZER} -fsanitize=undefined")

add_compile_options (${SANITIZER})
add_compile_options (-g -O3 -Wall -Wno-conversion)

include (FetchContent)

# re-extract downloaded archives when their URL changes
if (POLICY CMP0135)
  cmake_policy (SET CMP0135 NEW)
endif ()

#
# Fortran Expokit library and samples.  EXPOKIT is not bundled: it is
# fetched here, and it is non-commercial-use-only.  See NOTICE.
#
if (EXPOKIT_BUILD_FORTRAN)
  enable_language (Fortran)
  # only the Fortran library links them
  find_package (BLAS REQUIRED)
  find_package (LAPACK REQUIRED)

  FetchContent_Declare (
    expokit
    URL https://www.maths.uq.edu.au/expokit/expokit.tar.gz
    URL_HASH SHA256=aa9fc755f36349f3bdb9c9fc554bcb6184077621622e4d25c492751a3bfd8606
    DOWNLOAD_DIR "$ENV{HOME}/Downloads"
    SOURCE_SUBDIR "dummy"
    )
  FetchContent_MakeAvailable (expokit)
  FetchContent_GetProperties (expokit SOURCE_DIR EXPOKIT_SOURCE_DIR)

  # sample_d.f writes its FORMAT width with the <n> variable-format-expression
  # extension, which gfortran rejects; the file's own fixed-width FORMAT is
  # sitting right below it, commented out.  Swap the two.
  file (READ "${EXPOKIT_SOURCE_DIR}/fortran/sample_d.f" EXPOKIT_SAMPLE_D_SRC)
  string (REPLACE " 9001 format( <mprint>" "* 9001 format( <mprint>"
    EXPOKIT_SAMPLE_D_SRC "${EXPOKIT_SAMPLE_D_SRC}")
  string (REPLACE "* 9001 format( 5(1X" " 9001 format( 5(1X"
    EXPOKIT_SAMPLE_D_SRC "${EXPOKIT_SAMPLE_D_SRC}")
  set (EXPOKIT_SAMPLE_D "${CMAKE_CURRENT_BINARY_DIR}/expokit-patched/sample_d.f")
  file (WRITE "${EXPOKIT_SAMPLE_D}" "${EXPOKIT_SAMPLE_D_SRC}")

  add_library (expokitf
    ${EXPOKIT_SOURCE_DIR}/fortran/expokit.f
    ${EXPOKIT_SOURCE_DIR}/fortran/clock.f
    ${EXPOKIT_SOURCE_DIR}/fortran/mataid.f)
  target_link_libraries (expokitf ${LAPACK_LIBRARIES})
  target_compile_options (expokitf PRIVATE -fallow-argument-mismatch)

  foreach (EX g b z m d p)
    set (EXPOKIT_SAMPLE_SRC "${EXPOKIT_SOURCE_DIR}/fortran/sample_${EX}.f")
    if ("${EX}" STREQUAL "d")
      set (EXPOKIT_SAMPLE_SRC "${EXPOKIT_SAMPLE_D}")
    endif ()
    add_executable (sample_${EX} ${EXPOKIT_SAMPLE_SRC})
    #target_compile_options (sample_${EX} PRIVATE -fallow-argument-mismatch)
    target_link_libraries (sample_${EX} expokitf)
  endforeach (EX)
endif ()

if (NOT EXPOKIT_BUILD_TESTING)
  return ()
endif ()

#
# Tests
#
enable_testing ()

# Configure google-test as a downloadable library.
FetchContent_Declare (
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.18.0.tar.gz
  URL_HASH MD5=ae87634285dbe148a7abc6d19c4c05bd
  DOWNLOAD_DIR "$ENV{HOME}/Downloads"
  )

#
# Get Eigen
#
set (EXPOKIT_EIGEN_VERSION "5.0.1" CACHE STRING "Eigen version to build against")
set_property (CACHE EXPOKIT_EIGEN_VERSION PROPERTY STRINGS "3.3.7" "3.3.8" "5.0.1")

if (EXPOKIT_EIGEN_VERSION STREQUAL "3.3.7")
  set (EIGEN_MD5 b9e98a200d2455f06db9c661c5610496)
elseif (EXPOKIT_EIGEN_VERSION STREQUAL "3.3.8")
  set (EIGEN_MD5 432ef01499d514f4606343276afa0ec3)
elseif (EXPOKIT_EIGEN_VERSION STREQUAL "5.0.1")
  set (EIGEN_MD5 457323cbf3688f52b3623c0aa85044b7)
else ()
  message (FATAL_ERROR "Unknown EXPOKIT_EIGEN_VERSION '${EXPOKIT_EIGEN_VERSION}',"
                       " expected 3.3.7, 3.3.8 or 5.0.1")
endif ()
message (STATUS "Using Eigen ${EXPOKIT_EIGEN_VERSION}")

set (EIGEN_URL
  https://gitlab.com/libeigen/eigen/-/archive/${EXPOKIT_EIGEN_VERSION}/eigen-${EXPOKIT_EIGEN_VERSION}.tar.bz2)
FetchContent_Declare (
  libeigen
  URL ${EIGEN_URL}
  URL_HASH MD5=${EIGEN_MD5}
  DOWNLOAD_DIR "$ENV{HOME}/Downloads"
  SOURCE_SUBDIR "dummy"
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set (gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set (INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable (googletest libeigen)
include (GoogleTest)

# Get eigen header location
FetchContent_GetProperties (libeigen SOURCE_DIR EIGEN_SOURCE_DIR)
include_directories (SYSTEM "${EIGEN_SOURCE_DIR}")

#
# Example C++ Expokit programs, one per scalar type - instantiating them
# all in one translation unit needs more memory than a small CI runner has.
#
foreach (EX float cfloat double cdouble)
  add_executable (example_${EX} eexpokit/example_${EX}.cpp)
  target_link_libraries (example_${EX} eigen-expokit::eexpokit ${SANLIB})
endforeach (EX)

#
# Test binaries
#
# build handle for CI, which builds the tests without the examples
add_custom_target (expokit-tests)

option (EXPOKIT_BUILD_DUAL_TESTS "Build the AD-scalar tests, which fetch cppduals" ON)

set (ALL_TESTS
  eexpokit/test_expokithpp
  eexpokit/test_expv
  eexpokit/test_api
  eexpokit/test_expmv
  eexpokit/test_mexpv
  eexpokit/test_phiv
  eexpokit/test_accuracy
  eexpokit/test_sparse
  # one scalar type per translation unit, for the same reason the examples are
  # split that way: two of them together needs 2.2GB of compiler memory
  eexpokit/test_scalars_float
  eexpokit/test_scalars_double
  eexpokit/test_scalars_ldouble
  eexpokit/test_scalars_cfloat
  eexpokit/test_scalars_cdouble
  # one translation unit per public header, so a header that only compiles
  # after somebody else's include is a build failure
  eexpokit/test_hdr_chbv
  eexpokit/test_hdr_expokit
  eexpokit/test_hdr_expv
  eexpokit/test_hdr_expmv
  eexpokit/test_hdr_loadcrs
  eexpokit/test_hdr_mexpv
  eexpokit/test_hdr_padm
  eexpokit/test_hdr_phiv
  )

# loadcrs() has no other coverage: the examples generate their own matrix,
# and the only CRS file to hand is the one EXPOKIT ships.
if (EXPOKIT_BUILD_FORTRAN)
  list (APPEND ALL_TESTS eexpokit/test_loadcrs)
endif ()

# cppduals supplies the forward-mode AD scalar.  OFF drops that coverage, and
# the download with it, for offline or minimal builds.
if (EXPOKIT_BUILD_DUAL_TESTS)
  FetchContent_Declare (
    cppduals
    URL https://gitlab.com/tesch1/cppduals/-/archive/v0.10.0/cppduals-v0.10.0.tar.bz2
    URL_HASH MD5=348b9ccb8b05958d1c759ef30bc7434a
    DOWNLOAD_DIR "$ENV{HOME}/Downloads"
    SOURCE_SUBDIR "dummy"
    )
  FetchContent_MakeAvailable (cppduals)
  FetchContent_GetProperties (cppduals SOURCE_DIR CPPDUALS_SOURCE_DIR)
  list (APPEND ALL_TESTS eexpokit/test_scalars_dual)
endif ()

foreach (TEST_PATH ${ALL_TESTS})
  string (REPLACE "/" "_" TEST ${TEST_PATH})
  add_executable (${TEST} ${TEST_PATH}.cpp)
  target_link_libraries (${TEST} eigen-expokit::eexpokit gtest_main)
  add_dependencies (expokit-tests ${TEST})
  gtest_discover_tests (${TEST} TEST_LIST ${TEST}_targets)
  set_tests_properties (${${TEST}_targets} PROPERTIES TIMEOUT 10)
endforeach (TEST_PATH)

if (EXPOKIT_BUILD_DUAL_TESTS)
  target_include_directories (eexpokit_test_scalars_dual SYSTEM PRIVATE
    "${CPPDUALS_SOURCE_DIR}")
endif ()

if (EXPOKIT_BUILD_FORTRAN)
  target_compile_definitions (eexpokit_test_loadcrs PRIVATE
    EEXPOKIT_CRS_FILE="${EXPOKIT_SOURCE_DIR}/data/c1024.crs")
endif ()
