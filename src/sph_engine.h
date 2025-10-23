#pragma once

#include "particle.h"
#include "spatial_hash.h"
#include "kernels.h"
#include <memory>
#include <vector>
#include <functional>
#include <chrono>

namespace sph {

// SPH simulation parameters
struct SPHParameters {
    float rest_density = 1000.0f;        // Rest density (kg/m³)
    float gas_constant = 2000.0f;        // Gas constant for pressure
    float viscosity = 0.001f;            // Dynamic viscosity
    float smoothing_length = 0.02f;      // Smoothing length (h)
    float particle_mass = 0.001f;        // Mass of each particle
    float timestep = 0.001f;             // Fixed timestep (s)
    float gravity = -9.81f;              // Gravity acceleration
    float damping = 0.99f;               // Velocity damping
    float CFL_factor = 0.4f;             // CFL condition factor

    // Boundary conditions
    struct Boundaries {
        float xmin = -1.0f, xmax = 1.0f;
        float ymin = -1.0f, ymax = 1.0f;
        float zmin = -1.0f, zmax = 1.0f;
    } bounds;

    // Neighbor search
    float neighbor_search_radius = 0.04f;  // 2 * smoothing_length
};

// SPH Engine - Main simulation class
class SPHEngine {
private:
    // Core components
    ParticleSystem particles_;
    std::unique_ptr<SpatialHash> spatial_hash_;
    std::unique_ptr<Kernel> kernel_;
    SPHParameters params_;

    // Simulation state
    float current_time_ = 0.0f;
    size_t step_count_ = 0;
    bool initialized_ = false;

    // Performance monitoring
    struct PerformanceStats {
        double total_time = 0.0;
        double neighbor_search_time = 0.0;
        double density_computation_time = 0.0;
        double force_computation_time = 0.0;
        double integration_time = 0.0;
        size_t max_neighbors = 0;
        size_t total_neighbor_queries = 0;
    } perf_stats_;

    // Internal buffers for optimization
    std::vector<float> densities_;
    std::vector<float> pressures_;
    std::vector<glm::vec3> accelerations_;
    std::vector<std::vector<size_t>> neighbor_lists_;

    // Force computation methods
    void compute_densities();
    void compute_pressures();
    void compute_forces();
    void apply_external_forces();
    void apply_boundary_conditions();

    // Integration methods
    void integrate_euler(float dt);
    void integrate_verlet(float dt);
    void integrate_leapfrog(float dt);

    // CFL condition for adaptive timestep
    float compute_cfl_timestep() const;

    // Neighbor search optimization
    void update_neighbor_lists();

    // SIMD-optimized versions
    void compute_densities_simd();
    void compute_forces_simd();

public:
    // Constructor/Destructor
    explicit SPHEngine(size_t max_particles = 1000000);
    ~SPHEngine() = default;

    // Initialization
    void initialize(const SPHParameters& params);
    void initialize_dam_break();
    void initialize_fluid_drop();
    void initialize_granular_flow();

    // Particle management
    void add_particles(const std::vector<Particle>& particles);
    void clear_particles();

    // Simulation control
    void step(float dt = 0.0f);  // dt = 0 means use adaptive timestep
    void run_steps(size_t num_steps, bool adaptive_timestep = true);

    // Accessors
    const ParticleSystem& get_particles() const { return particles_; }
    const SPHParameters& get_parameters() const { return params_; }
    float get_current_time() const { return current_time_; }
    size_t get_step_count() const { return step_count_; }

    // Parameter modification
    void set_parameters(const SPHParameters& params);
    void set_gravity(float gravity) { params_.gravity = gravity; }
    void set_viscosity(float viscosity) { params_.viscosity = viscosity; }
    void set_smoothing_length(float h);

    // Boundary conditions
    void set_boundaries(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);

    // Performance monitoring
    const PerformanceStats& get_performance_stats() const { return perf_stats_; }
    void reset_performance_stats();

    // Utility functions
    void compute_conservation_errors(float& mass_error, float& energy_error) const;
    float get_total_mass() const;
    float get_total_energy() const;

    // Data export for visualization
    std::vector<glm::vec3> get_positions() const { return particles_.get_positions(); }
    std::vector<glm::vec3> get_velocities() const { return particles_.get_velocities(); }
    std::vector<float> get_densities() const { return densities_; }
    std::vector<float> get_pressures() const { return pressures_; }

    // Debug and validation
    void validate_simulation() const;
    bool is_initialized() const { return initialized_; }
};

// Integration methods
enum class IntegrationMethod {
    EULER,
    VERLET,
    LEAPFROG
};

// SPH equation implementations
namespace equations {

// Density computation: ρᵢ = ∑ⱼ mⱼ W(rᵢⱼ, h)
float compute_density(size_t i, const std::vector<size_t>& neighbors,
                     const ParticleSystem& particles, const Kernel& kernel);

// Pressure computation: Pᵢ = k(ρᵢ - ρ₀)
float compute_pressure(float density, float rest_density, float gas_constant);

// Pressure force: Fᵢ = -∑ⱼ mⱼ (Pᵢ + Pⱼ)/(2ρⱼ) ∇W(rᵢⱼ, h)
glm::vec3 compute_pressure_force(size_t i, const std::vector<size_t>& neighbors,
                                const ParticleSystem& particles, const Kernel& kernel,
                                const std::vector<float>& pressures, const std::vector<float>& densities);

// Viscosity force: Fᵢ = μ ∑ⱼ mⱼ (vⱼ - vᵢ)/ρⱼ ∇²W(rᵢⱼ, h)
glm::vec3 compute_viscosity_force(size_t i, const std::vector<size_t>& neighbors,
                                 const ParticleSystem& particles, const Kernel& kernel,
                                 float viscosity, const std::vector<float>& densities);

// External forces (gravity)
glm::vec3 compute_external_forces(float gravity);

// Total acceleration
glm::vec3 compute_acceleration(const glm::vec3& pressure_force,
                             const glm::vec3& viscosity_force,
                             const glm::vec3& external_force,
                             float particle_mass);

} // namespace equations

// Utility functions
namespace utils {

// Generate initial conditions for common scenarios
std::vector<Particle> create_dam_break_setup(const glm::vec3& dam_size,
                                           const glm::vec3& fluid_size,
                                           float spacing,
                                           const SPHParameters& params);

std::vector<Particle> create_fluid_drop_setup(const glm::vec3& center,
                                            float radius,
                                            float spacing,
                                            const SPHParameters& params);

std::vector<Particle> create_granular_flow_setup(const glm::vec3& pile_size,
                                               const glm::vec3& domain_size,
                                               float spacing,
                                               const SPHParameters& params);

// Validation functions
bool validate_particle_setup(const std::vector<Particle>& particles);
float compute_average_neighbors(const SPHEngine& engine);

} // namespace utils

} // namespace sph

