cmake_minimum_required(VERSION 3.22)

# Project metadata
project(CodyFortran LANGUAGES Fortran)
set(PROJECT_DESCRIPTION "CodyFortran: Correlated Dynamics with Fortran")
set(PROJECT_URL "https://github.com/fablackner/codyfortran")

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

# Optional cuBLAS backend (M_Utils_CublasLib, fieldFree.eigen.gpuQ):
# compile exactly one implementation submodule of M_Utils_CublasLib
option(CODY_USE_CUDA "Build the cuBLAS GPU backend (needs the CUDA toolkit, no CUDA compiler)" OFF)
if(CODY_USE_CUDA)
    list(FILTER sources EXCLUDE REGEX ".*/S_Utils_CublasLib_Stub\\.f90$")
else()
    list(FILTER sources EXCLUDE REGEX ".*/S_Utils_CublasLib_Cuda\\.f90$")
endif()

# Libxc exchange-correlation functionals (M_Utils_XcLib, sysXc: DFT and
# TDDFT): compile exactly one implementation submodule of M_Utils_XcLib
option(CODY_USE_LIBXC "Build with the Libxc exchange-correlation functionals (sysXc)" ON)
if(CODY_USE_LIBXC)
    list(FILTER sources EXCLUDE REGEX ".*/S_Utils_XcLib_Stub\\.f90$")
else()
    list(FILTER sources EXCLUDE REGEX ".*/S_Utils_XcLib_Libxc\\.f90$")
endif()

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

if(CODY_USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(CodyFortran PRIVATE CUDA::cublas CUDA::cudart)
endif()

if(CODY_USE_LIBXC)
    target_link_libraries(CodyFortran PRIVATE xc)
endif()

# Conditionally link compiler-specific libraries
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_link_libraries(CodyFortran PRIVATE fftw3 openblas)
endif()

# Umbrella targets for building all apps or all tests
add_custom_target(apps)
add_custom_target(tests)

# === Functions to compile Fortran programs ===

# One executable of the program `source`, written to `exe_dir` as <name>.exe.
# The target is named after the path of the executable (app/a/b/P_X ->
# app_a_b_P_X), so `cody build app/a/b/P_X` builds it.
function(add_program source exe_dir extra_libs umbrella_target)
    get_filename_component(exe_name "${source}" NAME_WE)
    file(RELATIVE_PATH rel_dir "${CMAKE_SOURCE_DIR}" "${exe_dir}")
    string(REPLACE "/" "_" target_name "${rel_dir}/${exe_name}")

    add_executable("${target_name}" "${source}")
    set_target_properties("${target_name}" PROPERTIES
        OUTPUT_NAME "${exe_name}.exe"
        RUNTIME_OUTPUT_DIRECTORY "${exe_dir}"
    )

    # Link CodyFortran and any extra libraries (like test-drive)
    target_link_libraries("${target_name}" PRIVATE CodyFortran ${extra_libs})

    # Add to umbrella target
    add_dependencies(${umbrella_target} "${target_name}")
endfunction()

# Every program source <name>.f90 below `dir` becomes an executable next to
# it. A program shared by several runs lies once in a common parent folder of
# their run folders: a folder that holds an input P_X.json but no P_X.f90
# gets P_X.exe from the nearest P_X.f90 above it (a program reads the JSON
# next to its executable), and such a shared source gets no executable of its
# own. Inputs in runs/ folders (generated, passed as argument) do not count.
function(create_executables_from_dir dir extra_libs umbrella_target)
    file(GLOB_RECURSE programs CONFIGURE_DEPENDS "${dir}/*.f90")
    file(GLOB_RECURSE inputs CONFIGURE_DEPENDS "${dir}/P_*.json")
    set(top_dir "${CMAKE_SOURCE_DIR}/${dir}")

    set(shared_programs "")
    foreach(input ${inputs})
        if(input MATCHES "/runs/")
            continue()
        endif()
        get_filename_component(run_dir "${input}" DIRECTORY)
        get_filename_component(name "${input}" NAME_WE)
        if(EXISTS "${run_dir}/${name}.f90")
            continue()
        endif()

        # Nearest program of that name in a parent folder
        set(source "")
        set(search_dir "${run_dir}")
        while(NOT source AND NOT search_dir STREQUAL top_dir)
            get_filename_component(search_dir "${search_dir}" DIRECTORY)
            if(EXISTS "${search_dir}/${name}.f90")
                set(source "${search_dir}/${name}.f90")
            endif()
        endwhile()
        if(NOT source)
            message(FATAL_ERROR "${input}: no program ${name}.f90 in this folder or a parent folder up to ${dir}/")
        endif()

        add_program("${source}" "${run_dir}" "${extra_libs}" ${umbrella_target})
        list(APPEND shared_programs "${source}")
    endforeach()

    foreach(file ${programs})
        if(file IN_LIST shared_programs)
            continue()
        endif()
        get_filename_component(src_dir "${file}" DIRECTORY)
        add_program("${file}" "${src_dir}" "${extra_libs}" ${umbrella_target})
    endforeach()
endfunction()

# Process programs in app/ (only link CodyFortran)
create_executables_from_dir(app "" apps)

# Process programs in test/ (link CodyFortran and test-drive)
create_executables_from_dir(test "test-drive" tests)

# === Define CTest tests ===
function(add_fortran_test file)
    get_filename_component(test_name "${file}" NAME_WE)
    
    # Reconstruct target name (same logic as create_executables_from_dir)
    file(RELATIVE_PATH rel_path "${CMAKE_SOURCE_DIR}" "${file}")
    string(REPLACE "/" "_" target_name "${rel_path}")
    string(REPLACE ".f90" "" target_name "${target_name}")

    # Get the source directory for working directory
    get_filename_component(src_dir "${file}" DIRECTORY)

    # A test that needs an optional library declares it in a line
    # "!> ctest-requires: libxc" of its source; without the library it is
    # not registered (it still compiles against the stub)
    file(STRINGS "${file}" requires_libxc REGEX "^!> ctest-requires: libxc")
    if(requires_libxc AND NOT CODY_USE_LIBXC)
        message(STATUS "Skipped test ${file} (needs CODY_USE_LIBXC)")
        return()
    endif()

    add_test(NAME "${test_name}"
             COMMAND "${target_name}"
             WORKING_DIRECTORY "${src_dir}")
    
    # Set environment variables for tests
    set_tests_properties("${test_name}" PROPERTIES
        ENVIRONMENT "CODY_PROJECT_DIR=${CMAKE_SOURCE_DIR}")

    # A test that must stop with a specific error declares the expected
    # output in a line "!> ctest-pass-regex: <regex>" of its source; it then
    # passes if and only if its output matches (the exit code is ignored)
    file(STRINGS "${file}" pass_regex REGEX "^!> ctest-pass-regex: ")
    if(pass_regex)
        string(REPLACE "!> ctest-pass-regex: " "" pass_regex "${pass_regex}")
        set_tests_properties("${test_name}" PROPERTIES PASS_REGULAR_EXPRESSION "${pass_regex}")
    endif()

    # Tests in the same folder may consume files written by earlier ones
    # (e.g. orb*.in, coeffs.in), so run them in sequence; folders stay parallel
    string(MAKE_C_IDENTIFIER "${src_dir}" dir_key)
    set(prev_test "${last_test_in_${dir_key}}")
    if(prev_test)
        set_tests_properties("${test_name}" PROPERTIES DEPENDS "${prev_test}")
    endif()
    set(last_test_in_${dir_key} "${test_name}" PARENT_SCOPE)
    message(STATUS "Added test ${file}")
endfunction()

file(GLOB_RECURSE test_sources CONFIGURE_DEPENDS "test/*.f90")
list(SORT test_sources)
foreach(test_file ${test_sources})
    add_fortran_test("${test_file}")
endforeach()
