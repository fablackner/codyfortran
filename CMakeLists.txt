cmake_minimum_required(VERSION 3.22)

# Project metadata
project(CodyFortran LANGUAGES Fortran)
set(PROJECT_DESCRIPTION "CodyFortran: Correlated Dynamics with Fortran")
set(PROJECT_URL "https://github.com/flackner/codyfortran")

# Enable testing early
include(CTest)
enable_testing()

# Set module directory to avoid cluttering the build dir
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Dynamically construct compiler-specific mod directory path
set(FORTRAN_COMPILER_ID "${CMAKE_Fortran_COMPILER_ID}")
set(FORTRAN_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
set(FORTRAN_MODDIR "${FORTRAN_COMPILER_ID}-${FORTRAN_COMPILER_VERSION}")
set(THIRD_PARTY_DIR "${CMAKE_SOURCE_DIR}/thirdParty/installation")
set(THIRD_PARTY_INCLUDE "${THIRD_PARTY_DIR}/include")

# Locate Fortran source files
# Added CONFIGURE_DEPENDS so CMake re-runs if files are added/removed
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "src/*.f90")

# Main library target
add_library(CodyFortran ${sources})

# Define include paths
set(CODY_INCLUDE_DIRS
    "${THIRD_PARTY_INCLUDE}"
    "${THIRD_PARTY_INCLUDE}/fortran_stdlib/${FORTRAN_MODDIR}/fortran_stdlib"
    "${THIRD_PARTY_INCLUDE}/fortran_stdlib/${FORTRAN_MODDIR}/blas"
    "${THIRD_PARTY_INCLUDE}/fortran_stdlib/${FORTRAN_MODDIR}/lapack"
    "${THIRD_PARTY_INCLUDE}/test-drive/${FORTRAN_MODDIR}"
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    list(APPEND CODY_INCLUDE_DIRS $ENV{MKLROOT}/include/fftw)
endif()

# Use target_include_directories for better encapsulation
# PUBLIC ensures that targets linking to CodyFortran also get these include paths
target_include_directories(CodyFortran PUBLIC ${CODY_INCLUDE_DIRS})

# Library search path
# Use target_link_directories instead of global link_directories
target_link_directories(CodyFortran PUBLIC "${THIRD_PARTY_DIR}/lib" "${THIRD_PARTY_DIR}/lib64")

# Link external libraries
target_link_libraries(CodyFortran PRIVATE
    arpack
    gsl
    expokit
    jsonfortran
    fortran_stdlib
)

# Conditionally link compiler-specific libraries
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_link_libraries(CodyFortran PRIVATE fftw3 openblas)
endif()

# === Function to compile Fortran programs ===
function(create_executables_from_dir dir extra_libs)
    file(GLOB_RECURSE programs CONFIGURE_DEPENDS "${dir}/*.f90")
    foreach(file ${programs})
        get_filename_component(exe_name "${file}" NAME_WE)
        get_filename_component(src_dir "${file}" DIRECTORY)
        
        add_executable("${exe_name}.exe" "${file}")
        # Link CodyFortran and any extra libraries (like test-drive)
        target_link_libraries("${exe_name}.exe" PRIVATE CodyFortran ${extra_libs})
        
        # Install to the source directory (preserving original behavior)
        install(TARGETS "${exe_name}.exe" DESTINATION "${src_dir}")
    endforeach()
endfunction()

# Process programs in app/ (only link CodyFortran)
create_executables_from_dir(app "")

# Process programs in test/ (link CodyFortran and test-drive)
create_executables_from_dir(test "test-drive")

# === Define CTest tests ===
function(add_fortran_test file)
    get_filename_component(test_name "${file}" NAME_WE)
    add_test(NAME "${test_name}"
             COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name}.exe")
endfunction()

file(GLOB_RECURSE test_sources CONFIGURE_DEPENDS "test/*.f90")
foreach(test_file ${test_sources})
    add_fortran_test("${test_file}")
    message(STATUS "Added test ${test_file}")
endforeach()
