cmake_minimum_required(VERSION 3.22)

project(codyLibs LANGUAGES Fortran CXX C)

# Set the installation path
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/installation" CACHE PATH "" FORCE)

# Include ExternalProject module
include(ExternalProject)

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/json-fortran")
    ExternalProject_Add(
        json-fortran
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/json-fortran"
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=$ENV{CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            -DSKIP_DOC_GEN=ON
            -DUSE_GNU_INSTALL_CONVENTION=ON
            -DCMAKE_Fortran_MODULE_DIRECTORY:PATH=include
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/stdlib")
    ExternalProject_Add(
        stdlib
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/stdlib"
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=$ENV{CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            -DBUILD_TESTING=OFF
            "-DCMAKE_Fortran_FLAGS=$ENV{FFLAGS} -DSTDLIB_EXTERNAL_BLAS -DSTDLIB_EXTERNAL_LAPACK"
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/test-drive")
    ExternalProject_Add(
        test-drive
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/test-drive"
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=$ENV{CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/FTorch")
    ExternalProject_Add(
        FTorch
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/FTorch"
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=$ENV{CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            -DBUILD_SHARED_LIBS=ON
            -DCMAKE_PREFIX_PATH=$ENV{HOME}/miniforge3/envs/torch/lib/python3.12/site-packages/torch
            -DTorch_DIR=$ENV{HOME}/miniforge3/envs/torch/lib/python3.12/site-packages/torch/share/cmake/Torch
            -DGPU_DEVICE=NONE
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/openblas")
    ExternalProject_Add(
        openblas
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/openblas"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -j4 PREFIX=${CMAKE_INSTALL_PREFIX} USE_OPENMP=1
        INSTALL_COMMAND make install PREFIX=${CMAKE_INSTALL_PREFIX} USE_OPENMP=1
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/fftw3")
    ExternalProject_Add(
        fftw3
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/fftw3"
        CONFIGURE_COMMAND
            <SOURCE_DIR>/configure
                --prefix=${CMAKE_INSTALL_PREFIX}
                --enable-openmp
                --enable-shared
                --enable-avx
        BUILD_COMMAND make -j4
        INSTALL_COMMAND make install
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/shtns")
    ExternalProject_Add(
        shtns
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/shtns"
        CONFIGURE_COMMAND
            <SOURCE_DIR>/configure
                --prefix=${CMAKE_INSTALL_PREFIX}
                --disable-python
                --enable-openmp
                --enable-shared
                --disable-f77
                $ENV{SHTNS_MKL_FLAG}
                --enable-march=native
                --enable-kernel-compiler=$ENV{CC}
        BUILD_COMMAND make -j4
        INSTALL_COMMAND make install PREFIX=${CMAKE_INSTALL_PREFIX} USE_OPENMP=1
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/expokit")
    ExternalProject_Add(
        expokit
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/expokit/fortran"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${CMAKE_Fortran_COMPILER} -fallow-argument-mismatch -c <SOURCE_DIR>/expokit.f -o expokit.o
        INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
            COMMAND ar rcs ${CMAKE_INSTALL_PREFIX}/lib/libexpokit.a expokit.o
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/arpack")
    ExternalProject_Add(
        arpack
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/arpack"
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=$ENV{CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
            -DMPI=OFF
            -DBUILD_SHARED_LIBS=ON
            -DICB=ON
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        BUILD_ALWAYS OFF
    )
endif()

#============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/sources/gsl")
    ExternalProject_Add(
        gsl
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/sources/gsl"
        CONFIGURE_COMMAND
            <SOURCE_DIR>/configure
                --prefix=${CMAKE_INSTALL_PREFIX}
        BUILD_COMMAND make -j4
        INSTALL_COMMAND make install
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS OFF
    )
endif()