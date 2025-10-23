#include "kernels.h"
#include <cmath>
#include <memory>
#include <stdexcept>
#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace sph {

// Base Kernel class
Kernel::Kernel(float h)
    : h_(h), h_inv_(1.0f / h), h_sq_(h * h) {
}

void Kernel::set_smoothing_length(float h) {
    h_ = h;
    h_inv_ = 1.0f / h;
    h_sq_ = h * h;
}

// Cubic Spline Kernel implementation
CubicSplineKernel::CubicSplineKernel(float h)
    : Kernel(h) {
    // Normalization constants for different dimensions
    // 3D: 1/(π h³)
    sigma_3d_ = 1.0f / (static_cast<float>(M_PI) * h * h * h);

    // 2D: 10/(7π h²)
    sigma_2d_ = 10.0f / (7.0f * static_cast<float>(M_PI) * h * h);

    // 1D: 2/(3h)
    sigma_1d_ = 2.0f / (3.0f * h);

    norm_ = sigma_3d_;
}

float CubicSplineKernel::W_1d(float q) const {
    if (q >= 0.0f && q <= 1.0f) {
        return sigma_1d_ * (2.0f/3.0f - q*q + 0.5f*q*q*q);
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        return sigma_1d_ * (1.0f/6.0f * tmp * tmp * tmp);
    }
    return 0.0f;
}

float CubicSplineKernel::W_2d(float q) const {
    if (q >= 0.0f && q <= 1.0f) {
        return sigma_2d_ * (2.0f/3.0f - q*q + 0.5f*q*q*q);
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        return sigma_2d_ * (1.0f/6.0f * tmp * tmp * tmp);
    }
    return 0.0f;
}

float CubicSplineKernel::W_3d(float q) const {
    if (q >= 0.0f && q <= 1.0f) {
        return sigma_3d_ * (2.0f/3.0f - q*q + 0.5f*q*q*q);
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        return sigma_3d_ * (1.0f/6.0f * tmp * tmp * tmp);
    }
    return 0.0f;
}

glm::vec3 CubicSplineKernel::gradW_1d(const glm::vec3& r, float q) const {
    float r_len = glm::length(r);
    if (r_len < 1e-6f) return glm::vec3(0.0f);

    glm::vec3 grad(0.0f);
    if (q >= 0.0f && q <= 1.0f) {
        grad = sigma_1d_ * (-2.0f * q + 1.5f * q * q) * (r / r_len);
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        grad = -sigma_1d_ * (0.5f * tmp * tmp) * (r / r_len);
    }
    return grad / h_;  // Chain rule: dW/dr = dW/dq * dq/dr = dW/dq * (r/(|r|h))
}

glm::vec3 CubicSplineKernel::gradW_2d(const glm::vec3& r, float q) const {
    float r_len = glm::length(r);
    if (r_len < 1e-6f) return glm::vec3(0.0f);

    glm::vec3 grad(0.0f);
    if (q >= 0.0f && q <= 1.0f) {
        grad = sigma_2d_ * (-2.0f * q + 1.5f * q * q) * (r / r_len);
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        grad = -sigma_2d_ * (0.5f * tmp * tmp) * (r / r_len);
    }
    return grad / h_;
}

glm::vec3 CubicSplineKernel::gradW_3d(const glm::vec3& r, float q) const {
    float r_len = glm::length(r);
    if (r_len < 1e-6f) return glm::vec3(0.0f);

    glm::vec3 grad(0.0f);
    if (q >= 0.0f && q <= 1.0f) {
        grad = sigma_3d_ * (-2.0f * q + 1.5f * q * q) * (r / r_len);
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        grad = -sigma_3d_ * (0.5f * tmp * tmp) * (r / r_len);
    }
    return grad / h_;
}

float CubicSplineKernel::laplacianW_1d(float q) const {
    if (q >= 0.0f && q <= 1.0f) {
        return sigma_1d_ * (-2.0f + 3.0f * q) / h_sq_;
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        return sigma_1d_ * tmp / h_sq_;
    }
    return 0.0f;
}

float CubicSplineKernel::laplacianW_2d(float q) const {
    if (q >= 0.0f && q <= 1.0f) {
        return sigma_2d_ * (-2.0f + 3.0f * q) / h_sq_;
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        return sigma_2d_ * tmp / h_sq_;
    }
    return 0.0f;
}

float CubicSplineKernel::laplacianW_3d(float q) const {
    if (q >= 0.0f && q <= 1.0f) {
        return sigma_3d_ * (-2.0f + 3.0f * q) / h_sq_;
    } else if (q > 1.0f && q <= 2.0f) {
        float tmp = 2.0f - q;
        return sigma_3d_ * tmp / h_sq_;
    }
    return 0.0f;
}

float CubicSplineKernel::W(const glm::vec3& r) const {
    float q = kernel_utils::compute_q(r, h_);
    return W_3d(q);
}

glm::vec3 CubicSplineKernel::gradW(const glm::vec3& r) const {
    float q = kernel_utils::compute_q(r, h_);
    return gradW_3d(r, q);
}

float CubicSplineKernel::laplacianW(const glm::vec3& r) const {
    float q = kernel_utils::compute_q(r, h_);
    return laplacianW_3d(q);
}

float CubicSplineKernel::W_2D(const glm::vec2& r) const {
    float q = glm::length(r) / h_;
    return W_2d(q);
}

float CubicSplineKernel::W_1D(float r) const {
    float q = std::abs(r) / h_;
    return W_1d(q);
}

// Wendland C2 Kernel implementation
WendlandC2Kernel::WendlandC2Kernel(float h)
    : Kernel(h) {
    norm_factor_ = 21.0f / (2.0f * static_cast<float>(M_PI) * h * h * h);
}

float WendlandC2Kernel::W(const glm::vec3& r) const {
    float q = kernel_utils::compute_q(r, h_);
    if (q >= 2.0f) return 0.0f;

    float tmp = 1.0f - 0.5f * q;
    return norm_factor_ * tmp * tmp * tmp * tmp * (2.0f * q + 1.0f);
}

glm::vec3 WendlandC2Kernel::gradW(const glm::vec3& r) const {
    float r_len = glm::length(r);
    if (r_len < 1e-6f) return glm::vec3(0.0f);

    float q = r_len / h_;
    if (q >= 2.0f) return glm::vec3(0.0f);

    float tmp = 1.0f - 0.5f * q;
    float dW_dq = -5.0f * tmp * tmp * tmp * q;

    return norm_factor_ * dW_dq * (r / (r_len * h_));
}

float WendlandC2Kernel::laplacianW(const glm::vec3& r) const {
    float q = kernel_utils::compute_q(r, h_);
    if (q >= 2.0f) return 0.0f;

    float tmp = 1.0f - 0.5f * q;
    return norm_factor_ * (5.0f / h_sq_) * tmp * tmp * (5.0f * q - 3.0f);
}

// Gaussian Kernel implementation
GaussianKernel::GaussianKernel(float h)
    : Kernel(h) {
    sigma_sq_inv_ = 1.0f / (h * h);
    norm_ = 1.0f / std::pow(static_cast<float>(M_PI) * h * h, 1.5f);
}

float GaussianKernel::W(const glm::vec3& r) const {
    float r_sq = glm::dot(r, r);
    return norm_ * std::exp(-r_sq * sigma_sq_inv_);
}

glm::vec3 GaussianKernel::gradW(const glm::vec3& r) const {
    float r_sq = glm::dot(r, r);
    float exp_term = std::exp(-r_sq * sigma_sq_inv_);
    return -2.0f * sigma_sq_inv_ * norm_ * exp_term * r;
}

float GaussianKernel::laplacianW(const glm::vec3& r) const {
    float r_sq = glm::dot(r, r);
    float exp_term = std::exp(-r_sq * sigma_sq_inv_);
    return 2.0f * sigma_sq_inv_ * norm_ * exp_term * (2.0f * sigma_sq_inv_ * r_sq - 3.0f);
}

// Kernel factory
std::unique_ptr<Kernel> create_kernel(KernelType type, float smoothing_length) {
    switch (type) {
        case KernelType::CUBIC_SPLINE:
            return std::make_unique<CubicSplineKernel>(smoothing_length);
        case KernelType::WENDLAND_C2:
            return std::make_unique<WendlandC2Kernel>(smoothing_length);
        case KernelType::GAUSSIAN:
            return std::make_unique<GaussianKernel>(smoothing_length);
        default:
            throw std::invalid_argument("Unknown kernel type");
    }
}

#ifdef __AVX2__

namespace avx {

// AVX2 cubic spline kernel for 8 simultaneous computations
__m256 cubic_spline_kernel_avx(__m256 distances, float h) {
    const __m256 zero = _mm256_setzero_ps();
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 two = _mm256_set1_ps(2.0f);

    // Compute q = distances / h
    __m256 h_vec = _mm256_set1_ps(h);
    __m256 q = _mm256_div_ps(distances, h_vec);

    // Normalization constant
    __m256 pi = _mm256_set1_ps(static_cast<float>(M_PI));
    __m256 h_cubed = _mm256_mul_ps(_mm256_mul_ps(h_vec, h_vec), h_vec);
    __m256 sigma = _mm256_div_ps(one, _mm256_mul_ps(pi, h_cubed));

    // q <= 1 branch
    __m256 mask_le1 = _mm256_cmp_ps(q, one, _CMP_LE_OQ);
    __m256 q_sq = _mm256_mul_ps(q, q);
    __m256 q_cu = _mm256_mul_ps(q_sq, q);
    __m256 w_le1 = _mm256_mul_ps(sigma, _mm256_add_ps(
        _mm256_sub_ps(_mm256_set1_ps(2.0f/3.0f), q_sq),
        _mm256_mul_ps(_mm256_set1_ps(0.5f), q_cu)
    ));

    // 1 < q <= 2 branch
    __m256 mask_gt1_le2 = _mm256_and_ps(
        _mm256_cmp_ps(q, one, _CMP_GT_OQ),
        _mm256_cmp_ps(q, two, _CMP_LE_OQ)
    );
    __m256 tmp = _mm256_sub_ps(two, q);
    __m256 tmp_cu = _mm256_mul_ps(_mm256_mul_ps(tmp, tmp), tmp);
    __m256 w_gt1_le2 = _mm256_mul_ps(sigma, _mm256_mul_ps(
        _mm256_set1_ps(1.0f/6.0f), tmp_cu
    ));

    // Combine results
    __m256 result = _mm256_blendv_ps(zero, w_le1, mask_le1);
    result = _mm256_blendv_ps(result, w_gt1_le2, mask_gt1_le2);

    return result;
}

__m256 cubic_spline_gradient_avx(__m256 distances, __m256 directions, float h) {
    // This is a simplified version - full implementation would be more complex
    // For now, return the kernel values (gradient implementation would need direction vectors)
    return cubic_spline_kernel_avx(distances, h);
}

} // namespace avx

#endif // __AVX2__

namespace kernel_utils {

float estimate_smoothing_length(float particle_spacing, float neighbor_count) {
    // Estimate based on desired neighbor count
    // For cubic spline kernel, typical neighbor count is ~50-100
    return particle_spacing * std::pow(neighbor_count / 50.0f, 1.0f/3.0f);
}

} // namespace kernel_utils

} // namespace sph

