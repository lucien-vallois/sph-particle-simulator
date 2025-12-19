#!/bin/bash

# SPH Particle Simulator - Cross-platform Build Script
# Author: Lucien Vallois
# Description: Automated build script for Linux/macOS/Windows (MSYS2)

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
BUILD_TYPE="Release"
BUILD_DIR="build"
USE_OPENMP="ON"
USE_CUDA="OFF"
BUILD_PYTHON="ON"
BUILD_EXAMPLES="ON"
BUILD_BENCHMARKS="ON"
BUILD_TESTS="OFF"
VERBOSE="OFF"
CLEAN_BUILD="OFF"

# Function to print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to detect OS
detect_os() {
    case "$(uname -s)" in
        Linux*)     OS="linux";;
        Darwin*)    OS="macos";;
        CYGWIN*|MINGW*|MSYS*) OS="windows";;
        *)          OS="unknown";;
    esac
}

# Function to check dependencies
check_dependencies() {
    print_info "Checking system dependencies..."

    # Check for required tools
    local missing_tools=()

    if ! command -v cmake &> /dev/null; then
        missing_tools+=("cmake")
    fi

    if ! command -v make &> /dev/null && ! command -v ninja &> /dev/null; then
        missing_tools+=("make/ninja")
    fi

    if [ "$OS" = "linux" ]; then
        if ! pkg-config --exists glfw3 2>/dev/null; then
            missing_tools+=("libglfw3-dev")
        fi
    elif [ "$OS" = "macos" ]; then
        if ! brew list glfw &> /dev/null; then
            missing_tools+=("glfw (brew install glfw)")
        fi
    fi

    if [ ${#missing_tools[@]} -ne 0 ]; then
        print_error "Missing dependencies: ${missing_tools[*]}"
        print_info "Please install them and run this script again."
        if [ "$OS" = "linux" ]; then
            print_info "Ubuntu/Debian: sudo apt install build-essential cmake libgl1-mesa-dev libglfw3-dev libglm-dev libomp-dev"
        elif [ "$OS" = "macos" ]; then
            print_info "macOS: brew install cmake glfw glm llvm"
        fi
        exit 1
    fi

    print_success "All dependencies found!"
}

# Function to parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--debug)
                BUILD_TYPE="Debug"
                shift
                ;;
            -r|--release)
                BUILD_TYPE="Release"
                shift
                ;;
            --no-openmp)
                USE_OPENMP="OFF"
                shift
                ;;
            --cuda)
                USE_CUDA="ON"
                shift
                ;;
            --no-python)
                BUILD_PYTHON="OFF"
                shift
                ;;
            --no-examples)
                BUILD_EXAMPLES="OFF"
                shift
                ;;
            --no-benchmarks)
                BUILD_BENCHMARKS="OFF"
                shift
                ;;
            --tests)
                BUILD_TESTS="ON"
                shift
                ;;
            -v|--verbose)
                VERBOSE="ON"
                shift
                ;;
            -c|--clean)
                CLEAN_BUILD="ON"
                shift
                ;;
            -h|--help)
                echo "SPH Particle Simulator Build Script"
                echo ""
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  -d, --debug          Build in Debug mode"
                echo "  -r, --release        Build in Release mode (default)"
                echo "  --no-openmp          Disable OpenMP parallelization"
                echo "  --cuda               Enable CUDA GPU acceleration"
                echo "  --no-python          Disable Python bindings"
                echo "  --no-examples        Don't build examples"
                echo "  --no-benchmarks      Don't build benchmarks"
                echo "  --tests              Build unit tests"
                echo "  -v, --verbose        Verbose output"
                echo "  -c, --clean          Clean build directory first"
                echo "  -h, --help           Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0                    # Standard release build"
                echo "  $0 --debug --tests    # Debug build with tests"
                echo "  $0 --cuda             # Release build with CUDA"
                echo "  $0 --clean            # Clean and rebuild"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                print_info "Use -h or --help for usage information"
                exit 1
                ;;
        esac
    done
}

# Function to clean build directory
clean_build() {
    if [ "$CLEAN_BUILD" = "ON" ] || [ ! -d "$BUILD_DIR" ]; then
        print_info "Cleaning build directory..."
        rm -rf "$BUILD_DIR"
        mkdir -p "$BUILD_DIR"
        print_success "Build directory cleaned"
    fi
}

# Function to configure CMake
configure_cmake() {
    print_info "Configuring CMake..."

    local cmake_args=(
        -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
        -DUSE_OPENMP="$USE_OPENMP"
        -DUSE_CUDA="$USE_CUDA"
        -DBUILD_PYTHON_BINDINGS="$BUILD_PYTHON"
        -DBUILD_EXAMPLES="$BUILD_EXAMPLES"
        -DBUILD_BENCHMARKS="$BUILD_BENCHMARKS"
        -DBUILD_TESTS="$BUILD_TESTS"
    )

    if [ "$VERBOSE" = "ON" ]; then
        cmake_args+=(-DCMAKE_VERBOSE_MAKEFILE=ON)
    fi

    cd "$BUILD_DIR"
    cmake .. "${cmake_args[@]}"

    if [ $? -eq 0 ]; then
        print_success "CMake configuration completed"
    else
        print_error "CMake configuration failed"
        exit 1
    fi

    cd ..
}

# Function to build project
build_project() {
    print_info "Building project..."

    cd "$BUILD_DIR"

    local make_cmd="make"
    if command -v ninja &> /dev/null; then
        make_cmd="ninja"
    fi

    local jobs=$(nproc 2>/dev/null || echo 4)

    if [ "$VERBOSE" = "ON" ]; then
        $make_cmd -j"$jobs" VERBOSE=1
    else
        $make_cmd -j"$jobs"
    fi

    if [ $? -eq 0 ]; then
        print_success "Build completed successfully!"
    else
        print_error "Build failed"
        exit 1
    fi

    cd ..
}

# Function to run tests
run_tests() {
    if [ "$BUILD_TESTS" = "ON" ]; then
        print_info "Running tests..."

        cd "$BUILD_DIR"

        if command -v ctest &> /dev/null; then
            ctest --output-on-failure
            if [ $? -eq 0 ]; then
                print_success "All tests passed!"
            else
                print_warning "Some tests failed"
            fi
        else
            print_warning "ctest not found, skipping tests"
        fi

        cd ..
    fi
}

# Function to run benchmarks
run_benchmarks() {
    if [ "$BUILD_BENCHMARKS" = "ON" ]; then
        print_info "Running benchmarks..."

        if [ -f "$BUILD_DIR/benchmarks/performance_test" ]; then
            "$BUILD_DIR/benchmarks/performance_test" 1000 5000 10000 25000
        else
            print_warning "Benchmark executable not found"
        fi
    fi
}

# Function to print build summary
print_summary() {
    print_success "Build Summary:"
    echo "  Build Type: $BUILD_TYPE"
    echo "  OpenMP: $USE_OPENMP"
    echo "  CUDA: $USE_CUDA"
    echo "  Python Bindings: $BUILD_PYTHON"
    echo "  Examples: $BUILD_EXAMPLES"
    echo "  Benchmarks: $BUILD_BENCHMARKS"
    echo "  Tests: $BUILD_TESTS"
    echo ""
    print_info "Executables created in: $BUILD_DIR/"
    if [ -f "$BUILD_DIR/sph_simulator" ]; then
        echo "  - sph_simulator (main executable)"
    fi
    if [ -d "$BUILD_DIR/examples" ]; then
        echo "  - examples/ (simulation examples)"
    fi
    if [ -d "$BUILD_DIR/benchmarks" ]; then
        echo "  - benchmarks/ (performance tests)"
    fi
    if [ "$BUILD_PYTHON" = "ON" ]; then
        echo "  - python/ (Python bindings)"
    fi
    echo ""
    print_info "To run the simulator:"
    echo "  cd $BUILD_DIR"
    echo "  ./sph_simulator"
    echo ""
    print_info "To run examples:"
    echo "  ./examples/dam_break 10000"
}

# Main execution
main() {
    detect_os
    print_info "SPH Particle Simulator Build Script"
    print_info "OS: $OS"
    echo ""

    parse_args "$@"
    check_dependencies
    clean_build
    configure_cmake
    build_project
    run_tests
    run_benchmarks
    print_summary

    print_success "Build process completed successfully! 🚀"
}

# Run main function with all arguments
main "$@"
