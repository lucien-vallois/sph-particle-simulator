#pragma once

#include <glm/glm.hpp>
#include <vector>

namespace sph {

// Smoothing kernel base class
class Kernel {
protected:
    float h_;        // Smoothing length
    float h_inv_;    // 1/h
    float h_sq_;     // h²
    float norm_;     // Normalization constant

public:
    explicit Kernel(float h);
    virtual ~Kernel() = default;

    // Core kernel functions
    virtual float W(const glm::vec3& r) const = 0;
    virtual glm::vec3 gradW(const glm::vec3& r) const = 0;
    virtual float laplacianW(const glm::vec3& r) const = 0;

    // Accessors
    float get_smoothing_length() const { return h_; }
    void set_smoothing_length(float h);

    // Utility functions
    float get_support_radius() const { return 2.0f * h_; }
    float get_cutoff_radius() const { return get_support_radius(); }
};

// Cubic spline kernel (most common in SPH)
class CubicSplineKernel : public Kernel {
private:
    float sigma_1d_;    // 1D normalization
    float sigma_2d_;    // 2D normalization
    float sigma_3d_;    // 3D normalization

    // Internal functions
    float W_1d(float q) const;
    float W_2d(float q) const;
    float W_3d(float q) const;

    glm::vec3 gradW_1d(const glm::vec3& r, float q) const;
    glm::vec3 gradW_2d(const glm::vec3& r, float q) const;
    glm::vec3 gradW_3d(const glm::vec3& r, float q) const;

    float laplacianW_1d(float q) const;
    float laplacianW_2d(float q) const;
    float laplacianW_3d(float q) const;

public:
    explicit CubicSplineKernel(float h);

    float W(const glm::vec3& r) const override;
    glm::vec3 gradW(const glm::vec3& r) const override;
    float laplacianW(const glm::vec3& r) const override;

    // Dimension-specific versions for special cases
    float W_3D(const glm::vec3& r) const { return W(r); }
    float W_2D(const glm::vec2& r) const;
    float W_1D(float r) const;
};

// Wendland C2 kernel (smoother, more stable)
class WendlandC2Kernel : public Kernel {
private:
    float norm_factor_;

public:
    explicit WendlandC2Kernel(float h);

    float W(const glm::vec3& r) const override;
    glm::vec3 gradW(const glm::vec3& r) const override;
    float laplacianW(const glm::vec3& r) const override;
};

// Gaussian kernel (for comparison, less common in SPH)
class GaussianKernel : public Kernel {
private:
    float sigma_sq_inv_;

public:
    explicit GaussianKernel(float h);

    float W(const glm::vec3& r) const override;
    glm::vec3 gradW(const glm::vec3& r) const override;
    float laplacianW(const glm::vec3& r) const override;
};

// Kernel factory
enum class KernelType {
    CUBIC_SPLINE,
    WENDLAND_C2,
    GAUSSIAN
};

std::unique_ptr<Kernel> create_kernel(KernelType type, float smoothing_length);

// SIMD-optimized kernel functions
#ifdef __AVX2__
namespace avx {

// AVX2 versions for processing 8 particles simultaneously
__m256 cubic_spline_kernel_avx(__m256 distances, float h);
__m256 cubic_spline_gradient_avx(__m256 distances, __m256 directions, float h);

} // namespace avx
#endif

// Utility functions for kernels
namespace kernel_utils {

// Compute q = |r| / h
inline float compute_q(const glm::vec3& r, float h) {
    return glm::length(r) / h;
}

// Check if position is within kernel support
inline bool within_support(const glm::vec3& r, float h, float support_factor = 2.0f) {
    return glm::dot(r, r) <= (support_factor * h) * (support_factor * h);
}

// Estimate optimal smoothing length based on particle spacing
float estimate_smoothing_length(float particle_spacing, float neighbor_count = 50.0f);

} // namespace kernel_utils

} // namespace sph

