#!/bin/bash

# Script: setBuildVars.sh
# Purpose: Sets environment variables for CMake-based Fortran build

# If you don't want pretty printed output from CODY, set this to OFF
export CODY_PRETTY_PRINT="ON"

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

# Set CMake build type (Capitalized for CMake convention)
if [[ "$BUILD_TYPE" == "debug" ]]; then
    export CMAKE_BUILD_TYPE="Debug"
else
    export CMAKE_BUILD_TYPE="Release"
fi

# Set Fortran compiler based on compiler family
if [[ "$COMPILER" == "gnu" ]]; then
    FC="gfortran"
else
    FC="ifx"
fi

# Set Fortran compiler
export FC=$FC

# Set compiler flags based on compiler and build type
GFORTRAN_DEBUG_FLAGS="-fopenmp -Werror -Wall -Wextra -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std \
-Wsurprising -Wunderflow -Wunused-parameter -Walign-commons -Wfunction-elimination -Wcharacter-truncation \
-Wline-truncation -Wtarget-lifetime -Wno-tabs -Winteger-division -pedantic -pedantic-errors -std=f2018 \
-fimplicit-none -fno-backslash -fall-intrinsics -frecursive -g -O0 -fbacktrace -fcheck=all \
-finit-real=snan -finit-integer=2147483647 -fmax-errors=10 -fstack-protector -fPIC -fprofile-arcs -ftest-coverage"

# maybe use -flto for link time optimization
# maybe use imprecise numerics -ffast-math
GFORTRAN_RELEASE_FLAGS="-fopenmp -O3 -march=native" 

# we should add -warn errors but then warning 5462 cannot be suppressed - fix in the future
INTEL_DEBUG_FLAGS="-qopenmp -warn all -debug all -check all -check nouninit -diag-enable all \
-standard-realloc-lhs -stand f18 -gen-interfaces -debug-parameters all -fstack-protector-all \
-g -traceback -O0 -qmkl"

# maybe use -flto for link time optimization
# maybe use -ipo for interprocedural optimization
# maybe use imprecise numerics -ffast-math -fp-model fast=2 -no-prec-div
INTEL_RELEASE_FLAGS="-qopenmp -O3 -xHost -qmkl"

# Set the appropriate flags based on compiler and build type
if [[ "$COMPILER" == "gnu" ]]; then
    if [[ "$BUILD_TYPE" = "debug" ]]; then
        FLAGS="${GFORTRAN_DEBUG_FLAGS}"
    elif [[ "$BUILD_TYPE" = "release" ]]; then
        FLAGS="${GFORTRAN_RELEASE_FLAGS}"
    else 
        echo "Error: Unrecognized build type. Use -t debug or -t release."
        return 1
    fi
elif [[ "$COMPILER" == "intel" ]]; then
    if [[ "$BUILD_TYPE" = "debug" ]]; then
        FLAGS="${INTEL_DEBUG_FLAGS}"
    elif [[ "$BUILD_TYPE" = "release" ]]; then
        FLAGS="${INTEL_RELEASE_FLAGS}"
    else 
        echo "Error: Unrecognized build type. Use -t debug or -t release."
        return 1
    fi
else
    echo "Error: Unrecognized compiler. Use -c gnu or -c intel."
    return 1
fi

# Set environment variables for flags imported by CMake
export FFLAGS="${FLAGS}"

# Set LD_LIBRARY_PATH to include thirdParty/lib for dynamic library resolution
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
THIRD_PARTY_LIB="${SCRIPT_DIR}/thirdParty/installation/lib:${SCRIPT_DIR}/thirdParty/installation/lib64"
export LD_LIBRARY_PATH="${THIRD_PARTY_LIB}:${LD_LIBRARY_PATH:-}"
export CODY_PROJECT_DIR="${SCRIPT_DIR}"

# Summary
echo "Build configuration set:"
echo "Compiler: $FC"
echo "Build Type: $BUILD_TYPE"
echo "Compiler Flags: ${FFLAGS}"
echo "Dynamic ld path: ${LD_LIBRARY_PATH}"
echo "CODY Project Dir: ${CODY_PROJECT_DIR}"
echo "CODY Pretty Print: ${CODY_PRETTY_PRINT}"