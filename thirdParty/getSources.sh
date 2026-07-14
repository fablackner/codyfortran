#!/bin/bash
set -e

# Set to "true" to download the library, "false" to skip
DOWNLOAD_EXPOKIT="true"
DOWNLOAD_SHTNS="true"
DOWNLOAD_STDLIB="true"
DOWNLOAD_TEST_DRIVE="true"
DOWNLOAD_JSON_FORTRAN="true"
DOWNLOAD_GSL="true"
DOWNLOAD_ARPACK="true"
# these are not needed when building with intel compilers and MKL
DOWNLOAD_FFTW3="true"
DOWNLOAD_OPENBLAS="true"
# these are not needed with out ml extensions
DOWNLOAD_FTORCH="false"

# Create sources and installation folders
mkdir -p sources
mkdir -p installation
mkdir -p installation/include

cd sources

# Clone or download required libraries
# expokit (Fortran version)
if [ "$DOWNLOAD_EXPOKIT" = "true" ] && [ ! -d "expokit" ]; then
    echo "Downloading expokit..."
    mkdir -p expokit
    curl -L -o expokit/expokit.tar.gz https://www.maths.uq.edu.au/expokit/expokit.tar.gz
    tar -zxvf expokit/expokit.tar.gz expokit/fortran/expokit.f
    rm -f expokit/expokit.tar.gz
fi

# FFTW3 (required by SHTNS with OpenMP support)
if [ "$DOWNLOAD_FFTW3" = "true" ] && [ ! -d "fftw3" ]; then
    echo "Downloading FFTW3..."
    curl -L -o fftw3.tar.gz https://www.fftw.org/fftw-3.3.10.tar.gz
    mkdir -p fftw3
    tar -xzf fftw3.tar.gz --strip-components=1 -C fftw3
    rm -f fftw3.tar.gz
fi

# SHTns
if [ "$DOWNLOAD_SHTNS" = "true" ] && [ ! -d "shtns" ]; then
    echo "Downloading SHTns..."
    git clone https://bitbucket.org/nschaeff/shtns.git
fi

# Fortran stdlib
if [ "$DOWNLOAD_STDLIB" = "true" ] && [ ! -d "stdlib" ]; then
    echo "Downloading Fortran stdlib..."
    git clone https://github.com/fortran-lang/stdlib.git
fi

# Fortran test-drive
if [ "$DOWNLOAD_TEST_DRIVE" = "true" ] && [ ! -d "test-drive" ]; then
    echo "Downloading Fortran test-drive..."
    git clone https://github.com/fortran-lang/test-drive.git
fi

# OpenBLAS
if [ "$DOWNLOAD_OPENBLAS" = "true" ] && [ ! -d "openblas" ]; then
    echo "Downloading OpenBLAS..."
    git clone https://github.com/xianyi/OpenBLAS.git openblas
fi

# ARPACK (ARnoldi PACKage for large-scale eigenvalue problems)
if [ "$DOWNLOAD_ARPACK" = "true" ] && [ ! -d "arpack" ]; then
    echo "Downloading ARPACK-NG..."
    git clone https://github.com/opencollab/arpack-ng.git arpack
fi

# JSON-Fortran
if [ "$DOWNLOAD_JSON_FORTRAN" = "true" ] && [ ! -d "json-fortran" ]; then
    echo "Downloading JSON-Fortran..."
    git clone https://github.com/jacobwilliams/json-fortran.git
fi

# GNU Scientific Library (GSL)
if [ "$DOWNLOAD_GSL" = "true" ] && [ ! -d "gsl" ]; then
    echo "Downloading GNU Scientific Library (GSL)..."
    curl -L -o gsl.tar.gz https://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
    mkdir -p gsl
    tar -xzf gsl.tar.gz --strip-components=1 -C gsl
    rm -f gsl.tar.gz
fi

# FTorch
if [ "$DOWNLOAD_FTORCH" = "true" ] && [ ! -d "FTorch" ]; then
    echo "Downloading FTorch..."
    git clone https://github.com/Cambridge-ICCS/FTorch.git
fi

echo "All requested libraries fetched into sources/"
