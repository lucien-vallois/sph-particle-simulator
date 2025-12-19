# SPH Particle Simulator - Docker Build
FROM ubuntu:22.04

# Labels for better maintainability
LABEL maintainer="Lucien Vallois <contact@example.com>"
LABEL description="High-performance SPH fluid simulation engine"
LABEL version="1.0.0"

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Build tools
    build-essential \
    cmake \
    ninja-build \
    git \
    # Graphics libraries
    libgl1-mesa-dev \
    libglfw3-dev \
    libglm-dev \
    libglew-dev \
    libglu1-mesa-dev \
    # OpenMP for parallelization
    libomp-dev \
    # X11 for GUI applications (if needed)
    libx11-dev \
    libxrandr-dev \
    libxinerama-dev \
    libxcursor-dev \
    libxi-dev \
    # Python for bindings and utilities
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-matplotlib \
    # Additional utilities
    wget \
    curl \
    vim \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip3 install --no-cache-dir \
    pybind11 \
    numpy \
    matplotlib \
    scipy \
    jupyter \
    notebook

# Set working directory
WORKDIR /app

# Copy source code
COPY . .

# Create build directory
RUN mkdir -p build

# Configure and build the project
RUN cd build && \
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DUSE_OPENMP=ON \
        -DBUILD_PYTHON_BINDINGS=ON \
        -DBUILD_EXAMPLES=ON \
        -DBUILD_BENCHMARKS=ON \
        -G Ninja && \
    ninja -j$(nproc)

# Create a non-root user for security
RUN useradd -m -s /bin/bash sphuser && \
    chown -R sphuser:sphuser /app
USER sphuser

# Set environment variables
ENV PYTHONPATH=/app/build/python:$PYTHONPATH
ENV DISPLAY=:0

# Create volume mount points
VOLUME ["/app/data", "/app/output"]

# Expose port for potential web interface (future use)
EXPOSE 8888

# Default command - run interactive shell
CMD ["/bin/bash"]

# Alternative commands for different use cases:
# CMD ["./build/sph_simulator"]                    # Run main simulator
# CMD ["./build/examples/dam_break", "10000"]      # Run dam break example
# CMD ["python3", "python/visualize.py"]           # Run Python visualization
# CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser"]  # Jupyter notebook
