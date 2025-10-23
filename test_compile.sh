#!/bin/bash

# Simple compilation test script for SPH simulator

echo "Testing SPH Simulator compilation..."

# Create build directory
mkdir -p build_test
cd build_test

# Configure with minimal options
cmake .. \
    -DBUILD_PYTHON_BINDINGS=OFF \
    -DUSE_CUDA=OFF \
    -DUSE_OPENMP=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_BENCHMARKS=OFF \
    -DCMAKE_BUILD_TYPE=Debug

# Try to compile
if make -j4; then
    echo "Compilation successful!"
    exit 0
else
    echo "Compilation failed!"
    exit 1
fi

