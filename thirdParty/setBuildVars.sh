#!/bin/bash

# Script: setBuildVars.sh
# Purpose: Sets environment variables for CMake-based Fortran build

# Detect if the script is being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: This script must be sourced, not executed."
    echo "Usage: source setBuildVars.sh [-c compiler] [-t type]"
    exit 1
fi

# Function to show help
show_help() {
    echo "Usage: source setBuildVars.sh [-c compiler] [-t type]"
    echo "Options:"
    echo "    -c compiler  Set the compiler family (gnu/intel, default: gnu)"
    echo "    -t type      Set build type (debug/release, default: release)"
    echo
    echo "Note: This script must be sourced to set environment variables"
    return 0
}

# Show help if requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    return 0
fi

# Default values
COMPILER="gnu"
BUILD_TYPE="release"

# Reset OPTIND because we might source this script multiple times
OPTIND=1

# Parse command line options
while getopts "c:t:h" opt; do
    case ${opt} in
        h)
            show_help
            return 0
            ;;
        c)
            COMPILER=$OPTARG
            ;;
        t)
            BUILD_TYPE=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            return 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument" 1>&2
            return 1
            ;;
    esac
done

# Validate compiler
if [[ "$COMPILER" != "gnu" && "$COMPILER" != "intel" ]]; then
    echo "Error: Invalid compiler family. Must be 'gnu' or 'intel'."
    return 1
fi

# Validate build type
if [[ "$BUILD_TYPE" != "debug" && "$BUILD_TYPE" != "release" ]]; then
    echo "Invalid build type. Must be 'debug' or 'release'."
    return 1
fi

# Set CMake build type
export CMAKE_BUILD_TYPE=$BUILD_TYPE

# Set compiler variables based on compiler family
if [[ "$COMPILER" == "gnu" ]]; then
    export FC="gfortran"
    export F77="gfortran"
    export CC="gcc"
    export CXX="g++"
    export SHTNS_MKL_FLAG=""
elif [[ "$COMPILER" == "intel" ]]; then
    export FC="ifx"
    export F77="ifx"
    export CC="icx"
    export CXX="icpx"
    export SHTNS_MKL_FLAG="--enable-mkl"
fi

# Get script directory for absolute paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
INSTALL_DIR="$SCRIPT_DIR/installation"
INCLUDE_FLAGS="-I$INSTALL_DIR/include"

# Define common flag patterns to reduce duplication
DEBUG_BASE="-g -O0"
RELEASE_BASE="-O3"

# Set compiler-specific flags
if [[ "$COMPILER" == "gnu" ]]; then
    OMP_FLAG="-fopenmp"
    if [[ "$BUILD_TYPE" == "debug" ]]; then
        F_SPECIFIC="-fbacktrace -fcheck=all"
        C_SPECIFIC="-Wall"
    else
        F_SPECIFIC="-march=native"
        C_SPECIFIC="-march=native"
    fi
elif [[ "$COMPILER" == "intel" ]]; then
    OMP_FLAG="-qopenmp"
    if [[ "$BUILD_TYPE" == "debug" ]]; then
        F_SPECIFIC="-traceback -check all"
        C_SPECIFIC="-Wall"
    else
        F_SPECIFIC="-xHost -qmkl"
        C_SPECIFIC="-xHost -qmkl"
    fi
fi

# Set build-type specific flags
if [[ "$BUILD_TYPE" == "debug" ]]; then
    BASE_FLAGS="$DEBUG_BASE"
else
    BASE_FLAGS="$RELEASE_BASE"
fi

# Set final flags
export FFLAGS="$OMP_FLAG $BASE_FLAGS $F_SPECIFIC $INCLUDE_FLAGS"
export CFLAGS="$OMP_FLAG $BASE_FLAGS $C_SPECIFIC $INCLUDE_FLAGS"

# Set FCFLAGS and CXXFLAGS to the same values as FFLAGS and CFLAGS
export FCFLAGS="$FFLAGS"
export CXXFLAGS="$CFLAGS"
export CPPFLAGS="$INCLUDE_FLAGS"
export LDFLAGS="-L$INSTALL_DIR/lib -L$INSTALL_DIR/lib64"

# Summary
echo "Build configuration set:"
echo "Compiler Family: $COMPILER"
echo "Fortran Compiler: $FC"
echo "C Compiler: $CC"
echo "C++ Compiler: $CXX"
echo "Build Type: $BUILD_TYPE"
echo "FFLAGS: $FFLAGS"
echo "FCFLAGS: $FCFLAGS"
echo "CFLAGS: $CFLAGS"
echo "CXXFLAGS: $CXXFLAGS"
echo "CPPFLAGS: $CPPFLAGS"
echo "LDFLAGS: $LDFLAGS"
