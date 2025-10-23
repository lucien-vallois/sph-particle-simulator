#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <memory>

namespace sph {

// Particle types for multi-phase simulations
enum class ParticleType {
    FLUID,
    BOUNDARY,
    SOLID
};

// Main particle structure
struct Particle {
    // Position and motion
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 acceleration;

    // SPH properties
    float density;
    float pressure;
    float mass;

    // Physical properties
    ParticleType type;
    float temperature;
    float viscosity;

    // For visualization and debugging
    glm::vec3 color;
    int id;

    // Constructor
    Particle()
        : position(0.0f), velocity(0.0f), acceleration(0.0f)
        , density(0.0f), pressure(0.0f), mass(1.0f)
        , type(ParticleType::FLUID), temperature(293.15f), viscosity(0.001f)
        , color(0.0f, 0.5f, 1.0f), id(-1) {}

    Particle(const glm::vec3& pos, float mass = 1.0f, ParticleType type = ParticleType::FLUID)
        : position(pos), velocity(0.0f), acceleration(0.0f)
        , density(0.0f), pressure(0.0f), mass(mass)
        , type(type), temperature(293.15f), viscosity(0.001f)
        , color(0.0f, 0.5f, 1.0f), id(-1) {}
};

// Container for particles with memory management
class ParticleSystem {
private:
    std::vector<Particle> particles_;
    size_t capacity_;

public:
    explicit ParticleSystem(size_t capacity = 1000000);
    ~ParticleSystem() = default;

    // Memory management
    void reserve(size_t capacity);
    void resize(size_t size);
    void clear();

    // Accessors
    Particle& operator[](size_t index) { return particles_[index]; }
    const Particle& operator[](size_t index) const { return particles_[index]; }
    size_t size() const { return particles_.size(); }
    size_t capacity() const { return capacity_; }
    bool empty() const { return particles_.empty(); }

    // Iterators
    auto begin() { return particles_.begin(); }
    auto end() { return particles_.end(); }
    auto begin() const { return particles_.cbegin(); }
    auto end() const { return particles_.cend(); }

    // Add particles
    void add_particle(const Particle& particle);
    void add_particles(const std::vector<Particle>& new_particles);

    // Remove particles
    void remove_particle(size_t index);
    void remove_particles(const std::vector<size_t>& indices);

    // Utility functions
    std::vector<glm::vec3> get_positions() const;
    std::vector<glm::vec3> get_velocities() const;
    std::vector<float> get_densities() const;
    std::vector<float> get_pressures() const;

    // Set properties for all particles
    void set_mass(float mass);
    void set_viscosity(float viscosity);
    void set_temperature(float temperature);

    // Boundary conditions
    void apply_boundary_conditions(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
};

// Utility functions
glm::vec3 generate_random_position(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
std::vector<Particle> create_fluid_block(const glm::vec3& center, const glm::vec3& size, float spacing, float mass = 1.0f);
std::vector<Particle> create_boundary_box(const glm::vec3& center, const glm::vec3& size, float spacing, float mass = 1.0f);

} // namespace sph

