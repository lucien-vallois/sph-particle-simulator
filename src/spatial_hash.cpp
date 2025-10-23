#include "spatial_hash.h"
#include <algorithm>
#include <cmath>
#include <limits>
#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace sph {

SpatialHash::SpatialHash(float cell_size)
    : cell_size_(cell_size), inv_cell_size_(1.0f / cell_size) {
}

void SpatialHash::build(const std::vector<glm::vec3>& positions) {
    clear();
    positions_cache_ = positions;

    for (size_t i = 0; i < positions.size(); ++i) {
        const auto& pos = positions[i];
        glm::ivec3 grid_coords = get_grid_coords(pos);
        uint64_t hash = hash_position(grid_coords.x, grid_coords.y, grid_coords.z);
        grid_[hash].push_back(i);
    }
}

std::vector<size_t> SpatialHash::query(const glm::vec3& position, float radius) const {
    return query_squared(position, radius * radius);
}

std::vector<size_t> SpatialHash::query_squared(const glm::vec3& position, float radius_squared) const {
    std::vector<size_t> neighbors;

    glm::ivec3 center_cell = get_grid_coords(position);
    int cell_radius = static_cast<int>(std::ceil(std::sqrt(radius_squared) * inv_cell_size_)) + 1;

    // Query neighboring cells
    for (int dx = -cell_radius; dx <= cell_radius; ++dx) {
        for (int dy = -cell_radius; dy <= cell_radius; ++dy) {
            for (int dz = -cell_radius; dz <= cell_radius; ++dz) {
                glm::ivec3 cell_coords = center_cell + glm::ivec3(dx, dy, dz);
                uint64_t hash = hash_position(cell_coords.x, cell_coords.y, cell_coords.z);

                auto it = grid_.find(hash);
                if (it != grid_.end()) {
                    for (size_t particle_idx : it->second) {
                        if (within_radius_squared(position, positions_cache_[particle_idx], radius_squared)) {
                            neighbors.push_back(particle_idx);
                        }
                    }
                }
            }
        }
    }

    return neighbors;
}

const std::vector<size_t>* SpatialHash::get_cell_particles(int x, int y, int z) const {
    uint64_t hash = hash_position(x, y, z);
    auto it = grid_.find(hash);
    return (it != grid_.end()) ? &it->second : nullptr;
}

void SpatialHash::clear() {
    grid_.clear();
    positions_cache_.clear();
}

size_t SpatialHash::get_max_particles_per_cell() const {
    size_t max_particles = 0;
    for (const auto& cell : grid_) {
        max_particles = std::max(max_particles, cell.second.size());
    }
    return max_particles;
}

#ifdef __AVX2__

void SpatialHashSIMD::query_multiple(const std::vector<glm::vec3>& query_positions,
                                   float radius,
                                   std::vector<std::vector<size_t>>& results) const {
    const float radius_squared = radius * radius;
    results.resize(query_positions.size());

    // Process 8 queries simultaneously with AVX2
    const size_t num_queries = query_positions.size();
    for (size_t base_idx = 0; base_idx < num_queries; base_idx += 8) {
        const size_t batch_size = std::min(size_t(8), num_queries - base_idx);

        // Load query positions
        __m256 query_x, query_y, query_z;
        if (batch_size >= 8) {
            query_x = _mm256_set_ps(
                query_positions[base_idx+7].x, query_positions[base_idx+6].x,
                query_positions[base_idx+5].x, query_positions[base_idx+4].x,
                query_positions[base_idx+3].x, query_positions[base_idx+2].x,
                query_positions[base_idx+1].x, query_positions[base_idx].x
            );
            query_y = _mm256_set_ps(
                query_positions[base_idx+7].y, query_positions[base_idx+6].y,
                query_positions[base_idx+5].y, query_positions[base_idx+4].y,
                query_positions[base_idx+3].y, query_positions[base_idx+2].y,
                query_positions[base_idx+1].y, query_positions[base_idx].y
            );
            query_z = _mm256_set_ps(
                query_positions[base_idx+7].z, query_positions[base_idx+6].z,
                query_positions[base_idx+5].z, query_positions[base_idx+4].z,
                query_positions[base_idx+3].z, query_positions[base_idx+2].z,
                query_positions[base_idx+1].z, query_positions[base_idx].z
            );
        }

        // For each particle in the system, check distance to all queries in batch
        for (size_t particle_idx = 0; particle_idx < positions_cache_.size(); ++particle_idx) {
            const glm::vec3& particle_pos = positions_cache_[particle_idx];

            __m256 particle_x = _mm256_set1_ps(particle_pos.x);
            __m256 particle_y = _mm256_set1_ps(particle_pos.y);
            __m256 particle_z = _mm256_set1_ps(particle_pos.z);

            // Compute squared distances
            __m256 dx = _mm256_sub_ps(query_x, particle_x);
            __m256 dy = _mm256_sub_ps(query_y, particle_y);
            __m256 dz = _mm256_sub_ps(query_z, particle_z);

            __m256 dx2 = _mm256_mul_ps(dx, dx);
            __m256 dy2 = _mm256_mul_ps(dy, dy);
            __m256 dz2 = _mm256_mul_ps(dz, dz);

            __m256 dist_squared = _mm256_add_ps(_mm256_add_ps(dx2, dy2), dz2);
            __m256 radius_squared_vec = _mm256_set1_ps(radius_squared);

            // Check which queries are within radius
            __m256 mask = _mm256_cmp_ps(dist_squared, radius_squared_vec, _CMP_LE_OQ);

            // Extract results
            int mask_bits = _mm256_movemask_ps(mask);
            for (int i = 0; i < batch_size; ++i) {
                if (mask_bits & (1 << i)) {
                    results[base_idx + i].push_back(particle_idx);
                }
            }
        }
    }
}

#endif // __AVX2__

namespace spatial_hash_utils {

float estimate_cell_size(float smoothing_length) {
    // Cell size should be approximately equal to smoothing length for optimal performance
    return smoothing_length;
}

std::vector<uint64_t> precompute_hashes(const std::vector<glm::vec3>& positions, float inv_cell_size) {
    std::vector<uint64_t> hashes;
    hashes.reserve(positions.size());

    for (const auto& pos : positions) {
        int x = static_cast<int>(std::floor(pos.x * inv_cell_size));
        int y = static_cast<int>(std::floor(pos.y * inv_cell_size));
        int z = static_cast<int>(std::floor(pos.z * inv_cell_size));

        uint64_t hash = 0;
        hash |= (static_cast<uint64_t>(x & 0x1FFFFF) << 42);
        hash |= (static_cast<uint64_t>(y & 0x1FFFFF) << 21);
        hash |= (static_cast<uint64_t>(z & 0x1FFFFF));
        hashes.push_back(hash);
    }

    return hashes;
}

} // namespace spatial_hash_utils

} // namespace sph

