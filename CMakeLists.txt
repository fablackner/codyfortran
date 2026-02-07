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
    "${THIRD_PARTY_INCLUDE}/fortran_stdlib/${FORTRAN_MODDIR}"
    "${THIRD_PARTY_INCLUDE}/test-drive/${FORTRAN_MODDIR}"
    #"${THIRD_PARTY_INCLUDE}/ftorch"
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
    fortran_stdlib_quadrature
    fortran_stdlib_sorting
    fortran_stdlib_linalg
    fortran_stdlib_stats
    fortran_stdlib_io
    fortran_stdlib_specialfunctions
    fortran_stdlib_bitsets
    fortran_stdlib_strings
    fortran_stdlib_intrinsics
    fortran_stdlib_math
    fortran_stdlib_core
    #ftorch
)

# Conditionally link compiler-specific libraries
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_link_libraries(CodyFortran PRIVATE fftw3 openblas)
endif()

# === Function to compile Fortran programs ===
function(create_executables_from_dir dir extra_libs)
    file(GLOB_RECURSE programs CONFIGURE_DEPENDS "${dir}/*.f90")
    foreach(file ${programs})
        # Generate unique target name based on relative path
        file(RELATIVE_PATH rel_path "${CMAKE_SOURCE_DIR}" "${file}")
        string(REPLACE "/" "_" target_name "${rel_path}")
        string(REPLACE ".f90" "" target_name "${target_name}")

        get_filename_component(exe_name "${file}" NAME_WE)
        get_filename_component(src_dir "${file}" DIRECTORY)
        
        # Calculate output directory to mirror source structure in build dir
        file(RELATIVE_PATH rel_dir "${CMAKE_SOURCE_DIR}" "${src_dir}")
        set(output_dir "${CMAKE_BINARY_DIR}/${rel_dir}")

        add_executable("${target_name}" "${file}")
        set_target_properties("${target_name}" PROPERTIES 
            OUTPUT_NAME "${exe_name}.exe"
            RUNTIME_OUTPUT_DIRECTORY "${output_dir}"
        )
        
        # Link CodyFortran and any extra libraries (like test-drive)
        target_link_libraries("${target_name}" PRIVATE CodyFortran ${extra_libs})
        
        # Install to the source directory (preserving original behavior)
        install(TARGETS "${target_name}" DESTINATION "${src_dir}")
    endforeach()
endfunction()

# Process programs in app/ (only link CodyFortran)
create_executables_from_dir(app "")

# Process programs in test/ (link CodyFortran and test-drive)
create_executables_from_dir(test "test-drive")

# === Define CTest tests ===
function(add_fortran_test file)
    get_filename_component(test_name "${file}" NAME_WE)
    
    # Reconstruct target name
    file(RELATIVE_PATH rel_path "${CMAKE_SOURCE_DIR}" "${file}")
    string(REPLACE "/" "_" target_name "${rel_path}")
    string(REPLACE ".f90" "" target_name "${target_name}")

    add_test(NAME "${test_name}"
             COMMAND "${target_name}")
endfunction()

file(GLOB_RECURSE test_sources CONFIGURE_DEPENDS "test/*.f90")
foreach(test_file ${test_sources})
    add_fortran_test("${test_file}")
    message(STATUS "Added test ${test_file}")
endforeach()
