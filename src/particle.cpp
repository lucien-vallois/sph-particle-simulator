#include "particle.h"
#include <algorithm>
#include <random>
#include <iostream>

namespace sph {

ParticleSystem::ParticleSystem(size_t capacity)
    : capacity_(capacity) {
    particles_.reserve(capacity_);
}

void ParticleSystem::reserve(size_t capacity) {
    capacity_ = capacity;
    particles_.reserve(capacity_);
}

void ParticleSystem::resize(size_t size) {
    particles_.resize(size);
}

void ParticleSystem::clear() {
    particles_.clear();
}

void ParticleSystem::add_particle(const Particle& particle) {
    if (particles_.size() < capacity_) {
        particles_.push_back(particle);
        particles_.back().id = static_cast<int>(particles_.size() - 1);
    } else {
        std::cerr << "Warning: Particle capacity exceeded\n";
    }
}

void ParticleSystem::add_particles(const std::vector<Particle>& new_particles) {
    for (const auto& particle : new_particles) {
        add_particle(particle);
    }
}

void ParticleSystem::remove_particle(size_t index) {
    if (index < particles_.size()) {
        particles_.erase(particles_.begin() + index);
        // Update IDs
        for (size_t i = index; i < particles_.size(); ++i) {
            particles_[i].id = static_cast<int>(i);
        }
    }
}

void ParticleSystem::remove_particles(const std::vector<size_t>& indices) {
    // Sort indices in descending order to avoid invalidating indices
    std::vector<size_t> sorted_indices = indices;
    std::sort(sorted_indices.rbegin(), sorted_indices.rend());

    for (size_t index : sorted_indices) {
        if (index < particles_.size()) {
            particles_.erase(particles_.begin() + index);
        }
    }

    // Update IDs
    for (size_t i = 0; i < particles_.size(); ++i) {
        particles_[i].id = static_cast<int>(i);
    }
}

std::vector<glm::vec3> ParticleSystem::get_positions() const {
    std::vector<glm::vec3> positions;
    positions.reserve(particles_.size());
    for (const auto& particle : particles_) {
        positions.push_back(particle.position);
    }
    return positions;
}

std::vector<glm::vec3> ParticleSystem::get_velocities() const {
    std::vector<glm::vec3> velocities;
    velocities.reserve(particles_.size());
    for (const auto& particle : particles_) {
        velocities.push_back(particle.velocity);
    }
    return velocities;
}

std::vector<float> ParticleSystem::get_densities() const {
    std::vector<float> densities;
    densities.reserve(particles_.size());
    for (const auto& particle : particles_) {
        densities.push_back(particle.density);
    }
    return densities;
}

std::vector<float> ParticleSystem::get_pressures() const {
    std::vector<float> pressures;
    pressures.reserve(particles_.size());
    for (const auto& particle : particles_) {
        pressures.push_back(particle.pressure);
    }
    return pressures;
}

void ParticleSystem::set_mass(float mass) {
    for (auto& particle : particles_) {
        particle.mass = mass;
    }
}

void ParticleSystem::set_viscosity(float viscosity) {
    for (auto& particle : particles_) {
        particle.viscosity = viscosity;
    }
}

void ParticleSystem::set_temperature(float temperature) {
    for (auto& particle : particles_) {
        particle.temperature = temperature;
    }
}

void ParticleSystem::apply_boundary_conditions(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    const float damping = 0.8f; // Energy damping factor

    for (auto& particle : particles_) {
        // X boundaries
        if (particle.position.x < xmin) {
            particle.position.x = xmin;
            particle.velocity.x *= -damping;
        } else if (particle.position.x > xmax) {
            particle.position.x = xmax;
            particle.velocity.x *= -damping;
        }

        // Y boundaries
        if (particle.position.y < ymin) {
            particle.position.y = ymin;
            particle.velocity.y *= -damping;
        } else if (particle.position.y > ymax) {
            particle.position.y = ymax;
            particle.velocity.y *= -damping;
        }

        // Z boundaries
        if (particle.position.z < zmin) {
            particle.position.z = zmin;
            particle.velocity.z *= -damping;
        } else if (particle.position.z > zmax) {
            particle.position.z = zmax;
            particle.velocity.z *= -damping;
        }
    }
}

// Utility functions
glm::vec3 generate_random_position(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<float> disx(xmin, xmax);
    std::uniform_real_distribution<float> disy(ymin, ymax);
    std::uniform_real_distribution<float> disz(zmin, zmax);

    return glm::vec3(disx(gen), disy(gen), disz(gen));
}

std::vector<Particle> create_fluid_block(const glm::vec3& center, const glm::vec3& size, float spacing, float mass) {
    std::vector<Particle> particles;

    int nx = static_cast<int>(size.x / spacing);
    int ny = static_cast<int>(size.y / spacing);
    int nz = static_cast<int>(size.z / spacing);

    glm::vec3 start = center - size * 0.5f;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                glm::vec3 position = start + glm::vec3(i * spacing, j * spacing, k * spacing);
                Particle particle(position, mass, ParticleType::FLUID);
                // Color based on position for visualization
                particle.color = glm::vec3(0.0f, 0.5f + 0.5f * (position.y - start.y) / size.y, 1.0f);
                particles.push_back(particle);
            }
        }
    }

    return particles;
}

std::vector<Particle> create_boundary_box(const glm::vec3& center, const glm::vec3& size, float spacing, float mass) {
    std::vector<Particle> particles;

    glm::vec3 half_size = size * 0.5f;
    glm::vec3 start = center - half_size;

    // Create 6 faces of the box
    auto add_face = [&](int axis1, int axis2, int axis3, float value, const glm::vec3& normal) {
        int n1 = static_cast<int>(size[axis1] / spacing);
        int n2 = static_cast<int>(size[axis2] / spacing);

        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                glm::vec3 position = start;
                position[axis1] = start[axis1] + i * spacing;
                position[axis2] = start[axis2] + j * spacing;
                position[axis3] = value;

                Particle particle(position, mass, ParticleType::BOUNDARY);
                particle.color = glm::vec3(0.5f, 0.5f, 0.5f); // Gray for boundaries
                particles.push_back(particle);
            }
        }
    };

    // Bottom face (z = start.z)
    add_face(0, 1, 2, start.z, glm::vec3(0, 0, -1));
    // Top face (z = start.z + size.z)
    add_face(0, 1, 2, start.z + size.z, glm::vec3(0, 0, 1));
    // Left face (x = start.x)
    add_face(1, 2, 0, start.x, glm::vec3(-1, 0, 0));
    // Right face (x = start.x + size.x)
    add_face(1, 2, 0, start.x + size.x, glm::vec3(1, 0, 0));
    // Front face (y = start.y)
    add_face(0, 2, 1, start.y, glm::vec3(0, -1, 0));
    // Back face (y = start.y + size.y)
    add_face(0, 2, 1, start.y + size.y, glm::vec3(0, 1, 0));

    return particles;
}

} // namespace sph

