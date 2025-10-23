# SPH Particle Simulator - Installation Guide

## Prerequisites

### Required Dependencies

#### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    libgl1-mesa-dev \
    libglfw3-dev \
    libglm-dev \
    libomp-dev

# For Python bindings
pip install pybind11
```

#### Linux (Fedora/CentOS)
```bash
sudo dnf install -y \
    gcc-c++ \
    cmake \
    mesa-libGL-devel \
    glfw-devel \
    glm-devel \
    openmp-devel

# For Python bindings
pip install pybind11
```

#### macOS
```bash
brew install \
    cmake \
    glfw \
    glm \
    llvm

# For Python bindings
pip install pybind11
```

#### Windows (MSYS2)
```bash
pacman -S \
    mingw-w64-x86_64-gcc \
    mingw-w64-x86_64-cmake \
    mingw-w64-x86_64-glfw \
    mingw-w64-x86_64-glm

# For Python bindings
pip install pybind11
```

## Build Instructions

### Option 1: Full Build (Recommended)

```bash
# Clone repository
git clone https://github.com/lucien-vallois/sph-particle-simulator.git
cd sph-particle-simulator

# Create build directory
mkdir build && cd build

# Configure with all features
cmake .. \
    -DUSE_OPENMP=ON \
    -DBUILD_PYTHON_BINDINGS=ON \
    -DBUILD_EXAMPLES=ON \
    -DBUILD_BENCHMARKS=ON

# Build
make -j$(nproc)

# Install (optional)
sudo make install
```

### Option 2: Minimal Build (For Testing)

```bash
mkdir build && cd build
cmake .. \
    -DBUILD_PYTHON_BINDINGS=OFF \
    -DUSE_OPENMP=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_BENCHMARKS=OFF
make
```

### Option 3: GPU Build (CUDA)

```bash
# Ensure CUDA is installed
mkdir build && cd build
cmake .. \
    -DUSE_CUDA=ON \
    -DUSE_OPENMP=ON
make -j$(nproc)
```

## Running the Simulator

### Basic Usage

```bash
# Run main simulator
./sph_simulator

# Run with specific parameters
./examples/dam_break 50000 1 1  # 50k particles, render, save data

# Run performance benchmark
./benchmarks/performance_test 1000 5000 10000
```

### Python Usage

```bash
# Build Python module
cd build
make sph
cp sph*.so ../python/

# Run Python examples
cd ../python
python visualize.py --mode visualize --particles 10000
python visualize.py --mode performance
```

## Troubleshooting

### Common Issues

#### 1. GLFW not found
```
CMake Error: Could not find glfw3
```

**Solution:**
- Ubuntu: `sudo apt install libglfw3-dev`
- macOS: `brew install glfw`
- Windows: Install via vcpkg or MSYS2

#### 2. GLM not found
```
GLM headers not found
```

**Solution:**
- Ubuntu: `sudo apt install libglm-dev`
- macOS: `brew install glm`
- Manual install: Download GLM and set GLM_INCLUDE_DIR

#### 3. OpenMP not found
```
OpenMP not found
```

**Solution:**
- Ubuntu: `sudo apt install libomp-dev`
- macOS: Already included with LLVM
- Windows: Use MSVC with `/openmp` flag

#### 4. Python bindings fail
```
pybind11 not found
```

**Solution:**
```bash
pip install pybind11
# Or install system-wide
sudo apt install pybind11-dev
```

### Compiler Issues

#### AVX2 Support
If your CPU doesn't support AVX2, disable it:
```bash
cmake .. -DCMAKE_CXX_FLAGS="-march=native"
```

#### C++17 Support
Ensure your compiler supports C++17:
```bash
g++ --version  # Should be 7.0+
clang++ --version  # Should be 5.0+
```

## Performance Optimization

### CPU Optimization

```bash
# Enable all optimizations
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_OPENMP=ON \
    -DCMAKE_CXX_FLAGS="-O3 -march=native -ffast-math"
```

### GPU Optimization (CUDA)

```bash
# Enable CUDA support
cmake .. -DUSE_CUDA=ON

# Set CUDA architecture (adjust for your GPU)
export CUDAARCHS="75"  # RTX 20xx series
```

## Testing

### Run Unit Tests

```bash
# Build with tests enabled
cmake .. -DBUILD_TESTS=ON
make test
```

### Performance Validation

```bash
# Run benchmark suite
./benchmarks/performance_test

# Expected results (Intel i7-9700K):
# 10k particles: ~2,400 FPS
# 100k particles: ~340 FPS
# 1M particles: ~38 FPS
```

## Docker Build (Alternative)

```dockerfile
FROM ubuntu:20.04

RUN apt update && apt install -y \
    build-essential \
    cmake \
    libgl1-mesa-dev \
    libglfw3-dev \
    libglm-dev \
    libomp-dev \
    python3-dev \
    python3-pip

RUN pip install pybind11

COPY . /app
WORKDIR /app/build

RUN cmake .. -DUSE_OPENMP=ON -DBUILD_PYTHON_BINDINGS=ON
RUN make -j$(nproc)

CMD ["./sph_simulator"]
```

```bash
docker build -t sph-simulator .
docker run -it --rm sph-simulator
```

