# 🤝 Contributing to SPH Particle Simulator

We welcome contributions to the SPH Particle Simulator! This document outlines the guidelines and best practices for contributing to this project.

## 📋 Table of Contents
- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)

## 🤝 Code of Conduct

This project follows a code of conduct to ensure a welcoming environment for all contributors. By participating, you agree to:

- Be respectful and inclusive
- Focus on constructive feedback
- Accept responsibility for mistakes
- Show empathy towards other contributors
- Help create a positive community

## 🚀 How to Contribute

### Types of Contributions

1. **🐛 Bug Reports**: Report bugs and issues
2. **💡 Feature Requests**: Suggest new features or improvements
3. **🔧 Code Contributions**: Submit pull requests with fixes or enhancements
4. **📚 Documentation**: Improve documentation, tutorials, or examples
5. **🧪 Testing**: Add or improve tests and benchmarks
6. **🎨 Design**: UI/UX improvements for visualization

### Getting Started

1. **Fork** the repository on GitHub
2. **Clone** your fork locally
3. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
4. **Make** your changes
5. **Test** thoroughly
6. **Commit** your changes (`git commit -m 'Add amazing feature'`)
7. **Push** to your branch (`git push origin feature/amazing-feature`)
8. **Open** a Pull Request

## 🛠️ Development Setup

### Prerequisites

**Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install -y build-essential cmake libgl1-mesa-dev libglfw3-dev libglm-dev libomp-dev python3-dev
pip install pybind11
```

**macOS:**
```bash
brew install cmake glfw glm llvm
pip install pybind11
```

**Windows:**
```powershell
# Using vcpkg or MSYS2
# Install: cmake, glfw3, glm, python3, pybind11
```

### Build Process

```bash
# Clone repository
git clone https://github.com/lucien-vallois/sph-particle-simulator.git
cd sph-particle-simulator

# Create build directory
mkdir build && cd build

# Configure (enable all features for development)
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DUSE_OPENMP=ON \
    -DBUILD_PYTHON_BINDINGS=ON \
    -DBUILD_EXAMPLES=ON \
    -DBUILD_BENCHMARKS=ON \
    -DBUILD_TESTS=ON

# Build
make -j$(nproc)

# Run tests
ctest --output-on-failure

# Install (optional)
sudo make install
```

## 💻 Coding Standards

### C++ Standards

- **Language**: C++17 standard
- **Formatting**: Follow the existing code style
- **Naming**: Use `snake_case` for variables/functions, `PascalCase` for classes/types
- **Documentation**: Use Doxygen-style comments for public APIs

```cpp
/**
 * @brief Compute density using SPH approximation
 * @param particles Particle system
 * @param smoothing_length Smoothing kernel radius
 * @return Computed density field
 */
std::vector<float> compute_density(const ParticleSystem& particles, float smoothing_length) {
    // Implementation here
}
```

### Python Standards

- Follow PEP 8 style guide
- Use type hints for function parameters and return values
- Include docstrings for all public functions

```python
def simulate_fluid_flow(particles: int, duration: float) -> np.ndarray:
    """
    Simulate fluid flow using SPH method.

    Args:
        particles: Number of particles to simulate
        duration: Simulation duration in seconds

    Returns:
        Array of particle positions over time
    """
    # Implementation here
    pass
```

### Commit Messages

Use clear, descriptive commit messages:

```
feat: add GPU acceleration for density computation
fix: resolve memory leak in spatial hashing
docs: update installation instructions for Windows
test: add benchmark for 1M particle simulation
```

## 🧪 Testing

### Running Tests

```bash
# Build with tests enabled
cmake .. -DBUILD_TESTS=ON
make
ctest --output-on-failure

# Run specific test
ctest -R test_sph_engine

# Run benchmarks
./benchmarks/performance_test 10000 50000 100000
```

### Adding Tests

- **Unit Tests**: Test individual components (kernels, spatial hash, etc.)
- **Integration Tests**: Test complete simulation workflows
- **Performance Tests**: Benchmark different optimization strategies
- **Regression Tests**: Ensure bug fixes don't break existing functionality

### Test Coverage

Aim for high test coverage, especially for:
- Core physics algorithms (SPH equations)
- Optimization routines (SIMD, OpenMP)
- Boundary conditions and stability
- Performance-critical code paths

## 📝 Submitting Changes

### Pull Request Process

1. **Update Documentation**: Ensure README and docs reflect your changes
2. **Add Tests**: Include tests for new features or bug fixes
3. **Update Changelog**: Document significant changes
4. **Self-Review**: Check your code before submitting
5. **Request Review**: Tag maintainers for review

### PR Template

```markdown
## Description
Brief description of the changes made.

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Performance improvement
- [ ] Documentation update
- [ ] Refactoring

## Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Performance benchmarks maintained
- [ ] Manual testing completed

## Screenshots (if applicable)
Add screenshots of visual changes.

## Checklist
- [ ] Code follows project style guidelines
- [ ] Documentation updated
- [ ] Tests added/updated
- [ ] All CI checks pass
```

## 🐛 Reporting Issues

### Bug Reports

Use the bug report template and include:

- **Description**: Clear description of the issue
- **Steps to Reproduce**: Detailed steps to reproduce the bug
- **Expected Behavior**: What should happen
- **Actual Behavior**: What actually happens
- **Environment**: OS, compiler, hardware specs
- **Logs/Error Messages**: Include relevant output

### Feature Requests

For feature requests, please:

- Check if the feature already exists or is planned
- Provide a clear use case and justification
- Consider implementation complexity
- Suggest a design approach if possible

## 🎯 Areas for Contribution

### High Priority
- GPU acceleration improvements (CUDA/HIP)
- Multi-threading optimizations
- Memory usage optimization
- Cross-platform compatibility

### Medium Priority
- New physical models (viscoelasticity, plasticity)
- Advanced boundary conditions
- Python API enhancements
- Documentation improvements

### Research Areas
- Multi-phase flows (fluid-solid coupling)
- Machine learning integration
- Adaptive mesh refinement
- Parallel computing (MPI)

## 📞 Getting Help

- **Issues**: Use GitHub Issues for bugs and feature requests
- **Discussions**: Use GitHub Discussions for questions and ideas
- **Email**: Contact maintainers directly for sensitive matters

## 🙏 Recognition

Contributors will be recognized in:
- Repository contributors list
- Changelog entries
- Project documentation
- Academic citations (if applicable)

Thank you for contributing to the SPH Particle Simulator! 🚀

