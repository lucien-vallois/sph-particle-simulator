#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <unordered_map>
#include <array>
#include <memory>

namespace sph {

// 3D spatial hash for efficient neighbor search
class SpatialHash {
private:
    float cell_size_;
    float inv_cell_size_;
    std::unordered_map<uint64_t, std::vector<size_t>> grid_;
    std::vector<glm::vec3> positions_cache_;

    // Hash function for 3D coordinates
    uint64_t hash_position(int x, int y, int z) const {
        // 64-bit hash: use 21 bits per coordinate (fits within ~2M range)
        uint64_t hash = 0;
        hash |= (static_cast<uint64_t>(x & 0x1FFFFF) << 42);
        hash |= (static_cast<uint64_t>(y & 0x1FFFFF) << 21);
        hash |= (static_cast<uint64_t>(z & 0x1FFFFF));
        return hash;
    }

    // Get grid coordinates for a position
    glm::ivec3 get_grid_coords(const glm::vec3& position) const {
        return glm::ivec3(
            static_cast<int>(std::floor(position.x * inv_cell_size_)),
            static_cast<int>(std::floor(position.y * inv_cell_size_)),
            static_cast<int>(std::floor(position.z * inv_cell_size_))
        );
    }

public:
    explicit SpatialHash(float cell_size = 0.02f);
    ~SpatialHash() = default;

    // Build the spatial hash from particle positions
    void build(const std::vector<glm::vec3>& positions);

    // Query neighbors within radius
    std::vector<size_t> query(const glm::vec3& position, float radius) const;

    // Query neighbors with squared radius (faster)
    std::vector<size_t> query_squared(const glm::vec3& position, float radius_squared) const;

    // Get all particles in a specific cell
    const std::vector<size_t>* get_cell_particles(int x, int y, int z) const;

    // Clear the hash
    void clear();

    // Get statistics
    size_t get_total_cells() const { return grid_.size(); }
    size_t get_max_particles_per_cell() const;

    // Cell size accessors
    float get_cell_size() const { return cell_size_; }
    void set_cell_size(float cell_size) {
        cell_size_ = cell_size;
        inv_cell_size_ = 1.0f / cell_size_;
    }

private:
    // Check if two positions are within radius squared
    static bool within_radius_squared(const glm::vec3& pos1, const glm::vec3& pos2, float radius_squared) {
        glm::vec3 diff = pos1 - pos2;
        return glm::dot(diff, diff) <= radius_squared;
    }

    // Get neighboring cells (27 cells for 3D)
    static std::array<glm::ivec3, 27> get_neighbor_offsets() {
        std::array<glm::ivec3, 27> offsets;
        int idx = 0;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    offsets[idx++] = glm::ivec3(dx, dy, dz);
                }
            }
        }
        return offsets;
    }
};

// Optimized version with SIMD support
#ifdef __AVX2__
class SpatialHashSIMD : public SpatialHash {
public:
    explicit SpatialHashSIMD(float cell_size = 0.02f) : SpatialHash(cell_size) {}

    // SIMD-optimized neighbor queries for multiple particles
    void query_multiple(const std::vector<glm::vec3>& query_positions,
                       float radius,
                       std::vector<std::vector<size_t>>& results) const;
};
#endif

// Utility functions for spatial hashing
namespace spatial_hash_utils {

// Estimate optimal cell size based on smoothing length
float estimate_cell_size(float smoothing_length);

// Precompute hash values for better performance
std::vector<uint64_t> precompute_hashes(const std::vector<glm::vec3>& positions, float inv_cell_size);

} // namespace spatial_hash_utils

} // namespace sph

