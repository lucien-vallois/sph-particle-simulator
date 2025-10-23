cmake_minimum_required(VERSION 3.16)
project(SPHSimulator VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(USE_OPENMP "Enable OpenMP parallelization" ON)
option(USE_CUDA "Enable CUDA GPU acceleration" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings with pybind11" ON)
option(BUILD_EXAMPLES "Build example applications" ON)
option(BUILD_BENCHMARKS "Build performance benchmarks" ON)

# Find packages
find_package(OpenGL REQUIRED)

# GLFW3 - try different methods
find_package(glfw3 QUIET)
if(NOT glfw3_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLFW3 REQUIRED glfw3)
    include_directories(${GLFW3_INCLUDE_DIRS})
    link_directories(${GLFW3_LIBRARY_DIRS})
endif()

# OpenGL extension loader (GLAD preferred, GLEW as fallback)
find_package(glad QUIET)
if(NOT glad_FOUND)
    find_package(GLEW QUIET)
    if(NOT GLEW_FOUND)
        message(WARNING "Neither GLAD nor GLEW found. OpenGL extension loading may fail.")
        message("Please install GLAD (recommended) or GLEW:")
        message("  Ubuntu: sudo apt install libglad-dev or libglew-dev")
        message("  Or download GLAD from: https://glad.dav1d.de/")
    endif()
endif()

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

if(USE_CUDA)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
endif()

if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11 REQUIRED)
endif()

# Include directories
include_directories(src)

# GLM - try to find it, if not found assume it's available system-wide
find_package(glm QUIET)
if(NOT glm_FOUND)
    message(STATUS "GLM not found via find_package, assuming system installation")
    # Check if GLM headers are available
    find_path(GLM_INCLUDE_DIR glm/glm.hpp)
    if(GLM_INCLUDE_DIR)
        include_directories(${GLM_INCLUDE_DIR})
        message(STATUS "Found GLM headers at: ${GLM_INCLUDE_DIR}")
    else()
        message(WARNING "GLM headers not found. Please install GLM library.")
    endif()
endif()

# Source files
set(SPH_SOURCES
    src/sph_engine.cpp
    src/spatial_hash.cpp
    src/kernels.cpp
    src/renderer.cpp
    src/particle.cpp
)

set(SPH_HEADERS
    src/sph_engine.h
    src/spatial_hash.h
    src/kernels.h
    src/renderer.h
    src/particle.h
)

# Main library
add_library(sph_core STATIC ${SPH_SOURCES} ${SPH_HEADERS})

target_link_libraries(sph_core PUBLIC
    OpenGL::GL
)

if(glfw3_FOUND)
    target_link_libraries(sph_core PUBLIC glfw)
else()
    target_link_libraries(sph_core PUBLIC ${GLFW3_LIBRARIES})
endif()

# Link OpenGL extension loader
if(glad_FOUND)
    target_link_libraries(sph_core PUBLIC glad::glad)
elseif(GLEW_FOUND)
    target_link_libraries(sph_core PUBLIC GLEW::GLEW)
endif()

if(USE_OPENMP)
    target_link_libraries(sph_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(sph_core PUBLIC USE_OPENMP)
endif()

target_compile_definitions(sph_core PUBLIC
    GLM_FORCE_RADIANS
    GLM_FORCE_DEPTH_ZERO_TO_ONE
)

# Main executable
add_executable(sph_simulator src/main.cpp)
target_link_libraries(sph_simulator PRIVATE sph_core)

# Python bindings
if(BUILD_PYTHON_BINDINGS)
    pybind11_add_module(sph python/bindings.cpp)
    target_link_libraries(sph PRIVATE sph_core)
endif()

# Examples
if(BUILD_EXAMPLES)
    add_executable(dam_break examples/dam_break.cpp)
    target_link_libraries(dam_break PRIVATE sph_core)

    add_executable(fluid_drop examples/fluid_drop.cpp)
    target_link_libraries(fluid_drop PRIVATE sph_core)

    add_executable(granular_flow examples/granular_flow.cpp)
    target_link_libraries(granular_flow PRIVATE sph_core)
endif()

# Benchmarks
if(BUILD_BENCHMARKS)
    add_executable(performance_test benchmarks/performance_test.cpp)
    target_link_libraries(performance_test PRIVATE sph_core)
endif()

# Compiler optimizations
if(MSVC)
    target_compile_options(sph_core PRIVATE /W4 /arch:AVX2 /fp:fast)
    target_compile_options(sph_simulator PRIVATE /W4 /arch:AVX2 /fp:fast)
else()
    target_compile_options(sph_core PRIVATE -Wall -Wextra -mavx2 -ffast-math -fopenmp)
    target_compile_options(sph_simulator PRIVATE -Wall -Wextra -mavx2 -ffast-math -fopenmp)
endif()

# Install targets
install(TARGETS sph_simulator
    RUNTIME DESTINATION bin
)

if(BUILD_PYTHON_BINDINGS)
    install(TARGETS sph
        LIBRARY DESTINATION lib/python
    )
endif()

# Copy shaders to build directory
configure_file(shaders/particle.vert ${CMAKE_BINARY_DIR}/shaders/particle.vert COPYONLY)
configure_file(shaders/particle.frag ${CMAKE_BINARY_DIR}/shaders/particle.frag COPYONLY)
