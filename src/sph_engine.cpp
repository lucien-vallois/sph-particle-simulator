#include "sph_engine.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace sph {

// SPHEngine implementation
SPHEngine::SPHEngine(size_t max_particles)
    : particles_(max_particles) {
    // Initialize with default parameters
    spatial_hash_ = std::make_unique<SpatialHash>(params_.smoothing_length);
    kernel_ = std::make_unique<CubicSplineKernel>(params_.smoothing_length);
}

void SPHEngine::initialize(const SPHParameters& params) {
    params_ = params;
    spatial_hash_->set_cell_size(params.neighbor_search_radius);
    kernel_ = create_kernel(KernelType::CUBIC_SPLINE, params.smoothing_length);

    // Resize internal buffers
    densities_.resize(particles_.capacity());
    pressures_.resize(particles_.capacity());
    accelerations_.resize(particles_.capacity());
    neighbor_lists_.resize(particles_.capacity());

    initialized_ = true;
    reset_performance_stats();
}

void SPHEngine::initialize_dam_break() {
    if (!initialized_) {
        initialize(SPHParameters{});
    }

    auto particles = utils::create_dam_break_setup(
        glm::vec3(0.4f, 0.6f, 0.8f),  // dam size
        glm::vec3(0.2f, 0.4f, 0.8f),  // fluid size
        0.01f,                         // spacing
        params_
    );

    particles_.clear();
    particles_.add_particles(particles);
}

void SPHEngine::initialize_fluid_drop() {
    if (!initialized_) {
        initialize(SPHParameters{});
    }

    auto particles = utils::create_fluid_drop_setup(
        glm::vec3(0.0f, 0.5f, 0.0f),  // center
        0.1f,                          // radius
        0.008f,                        // spacing
        params_
    );

    particles_.clear();
    particles_.add_particles(particles);
}

void SPHEngine::initialize_granular_flow() {
    if (!initialized_) {
        initialize(SPHParameters{});
    }

    auto particles = utils::create_granular_flow_setup(
        glm::vec3(0.3f, 0.4f, 0.8f),  // pile size
        glm::vec3(1.0f, 1.0f, 1.0f),  // domain size
        0.012f,                        // spacing
        params_
    );

    particles_.clear();
    particles_.add_particles(particles);
}

void SPHEngine::add_particles(const std::vector<Particle>& particles) {
    particles_.add_particles(particles);
}

void SPHEngine::clear_particles() {
    particles_.clear();
    current_time_ = 0.0f;
    step_count_ = 0;
}

void SPHEngine::step(float dt) {
    if (!initialized_ || particles_.size() == 0) return;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Use adaptive timestep if dt == 0
    if (dt <= 0.0f) {
        dt = compute_cfl_timestep();
    }

    // Neighbor search
    auto neighbor_start = std::chrono::high_resolution_clock::now();
    update_neighbor_lists();
    auto neighbor_end = std::chrono::high_resolution_clock::now();
    perf_stats_.neighbor_search_time +=
        std::chrono::duration<double>(neighbor_end - neighbor_start).count();

    // Density computation
    auto density_start = std::chrono::high_resolution_clock::now();
    compute_densities();
    auto density_end = std::chrono::high_resolution_clock::now();
    perf_stats_.density_computation_time +=
        std::chrono::duration<double>(density_end - density_start).count();

    // Pressure computation
    compute_pressures();

    // Force computation
    auto force_start = std::chrono::high_resolution_clock::now();
    compute_forces();
    auto force_end = std::chrono::high_resolution_clock::now();
    perf_stats_.force_computation_time +=
        std::chrono::duration<double>(force_end - force_start).count();

    // Integration
    auto integration_start = std::chrono::high_resolution_clock::now();
    integrate_leapfrog(dt);
    auto integration_end = std::chrono::high_resolution_clock::now();
    perf_stats_.integration_time +=
        std::chrono::duration<double>(integration_end - integration_start).count();

    // Boundary conditions
    apply_boundary_conditions();

    // Update time and counters
    current_time_ += dt;
    step_count_++;

    auto end_time = std::chrono::high_resolution_clock::now();
    perf_stats_.total_time +=
        std::chrono::duration<double>(end_time - start_time).count();
}

void SPHEngine::run_steps(size_t num_steps, bool adaptive_timestep) {
    for (size_t i = 0; i < num_steps; ++i) {
        step(adaptive_timestep ? 0.0f : params_.timestep);
    }
}

void SPHEngine::set_parameters(const SPHParameters& params) {
    params_ = params;
    spatial_hash_->set_cell_size(params.neighbor_search_radius);
    set_smoothing_length(params.smoothing_length);
}

void SPHEngine::set_smoothing_length(float h) {
    params_.smoothing_length = h;
    params_.neighbor_search_radius = 2.0f * h;
    spatial_hash_->set_cell_size(params_.neighbor_search_radius);
    kernel_ = create_kernel(KernelType::CUBIC_SPLINE, h);
}

void SPHEngine::set_boundaries(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    params_.bounds.xmin = xmin;
    params_.bounds.xmax = xmax;
    params_.bounds.ymin = ymin;
    params_.bounds.ymax = ymax;
    params_.bounds.zmin = zmin;
    params_.bounds.zmax = zmax;
}

void SPHEngine::reset_performance_stats() {
    perf_stats_ = PerformanceStats{};
}

void SPHEngine::compute_conservation_errors(float& mass_error, float& energy_error) const {
    float total_mass = get_total_mass();
    float initial_mass = particles_.size() * params_.particle_mass;
    mass_error = std::abs(total_mass - initial_mass) / initial_mass;

    // Energy conservation (simplified - would need initial energy for proper calculation)
    energy_error = 0.0f;  // TODO: Implement proper energy conservation check
}

float SPHEngine::get_total_mass() const {
    return std::accumulate(densities_.begin(), densities_.begin() + particles_.size(), 0.0f) *
           (params_.smoothing_length * params_.smoothing_length * params_.smoothing_length);
}

float SPHEngine::get_total_energy() const {
    float kinetic_energy = 0.0f;
    for (size_t i = 0; i < particles_.size(); ++i) {
        const auto& particle = particles_[i];
        float speed_sq = glm::dot(particle.velocity, particle.velocity);
        kinetic_energy += 0.5f * particle.mass * speed_sq;
    }
    return kinetic_energy;
}

// Private implementation methods
void SPHEngine::compute_densities() {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        densities_[i] = equations::compute_density(i, neighbor_lists_[i],
                                                 particles_, *kernel_);
        particles_[i].density = densities_[i];
    }
}

void SPHEngine::compute_pressures() {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        pressures_[i] = equations::compute_pressure(densities_[i],
                                                  params_.rest_density,
                                                  params_.gas_constant);
        particles_[i].pressure = pressures_[i];
    }
}

void SPHEngine::compute_forces() {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        glm::vec3 pressure_force = equations::compute_pressure_force(i, neighbor_lists_[i],
                                                                   particles_, *kernel_,
                                                                   pressures_, densities_);

        glm::vec3 viscosity_force = equations::compute_viscosity_force(i, neighbor_lists_[i],
                                                                     particles_, *kernel_,
                                                                     params_.viscosity, densities_);

        glm::vec3 external_force = equations::compute_external_forces(params_.gravity);

        accelerations_[i] = equations::compute_acceleration(pressure_force, viscosity_force,
                                                          external_force, particles_[i].mass);
    }
}

void SPHEngine::apply_external_forces() {
    // Already included in compute_forces()
}

void SPHEngine::apply_boundary_conditions() {
    particles_.apply_boundary_conditions(
        params_.bounds.xmin, params_.bounds.xmax,
        params_.bounds.ymin, params_.bounds.ymax,
        params_.bounds.zmin, params_.bounds.zmax
    );
}

void SPHEngine::integrate_euler(float dt) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        auto& particle = particles_[i];
        particle.velocity += accelerations_[i] * dt;
        particle.position += particle.velocity * dt;

        // Apply damping
        particle.velocity *= params_.damping;
    }
}

void SPHEngine::integrate_verlet(float dt) {
    // Verlet integration (position-based)
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        auto& particle = particles_[i];
        glm::vec3 new_position = 2.0f * particle.position - particle.position +
                                accelerations_[i] * dt * dt;

        particle.velocity = (new_position - particle.position) / dt;
        particle.position = new_position;

        // Apply damping
        particle.velocity *= params_.damping;
    }
}

void SPHEngine::integrate_leapfrog(float dt) {
    // Leapfrog integration (better energy conservation)
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        auto& particle = particles_[i];

        // Half velocity update
        particle.velocity += 0.5f * accelerations_[i] * dt;

        // Position update
        particle.position += particle.velocity * dt;

        // Half velocity update (will be completed in next step with new acceleration)
        particle.velocity += 0.5f * accelerations_[i] * dt;

        // Apply damping
        particle.velocity *= params_.damping;
    }
}

float SPHEngine::compute_cfl_timestep() const {
    float max_velocity = 0.0f;
    float max_force = 0.0f;

    for (size_t i = 0; i < particles_.size(); ++i) {
        float velocity = glm::length(particles_[i].velocity);
        max_velocity = std::max(max_velocity, velocity);

        float force = glm::length(accelerations_[i]) * particles_[i].mass;
        max_force = std::max(max_force, force);
    }

    // CFL condition based on velocity
    float dt_cfl = params_.CFL_factor * params_.smoothing_length /
                   (max_velocity + 1e-6f);  // Avoid division by zero

    // Force-based condition
    float dt_force = params_.CFL_factor * std::sqrt(params_.smoothing_length /
                                                    (glm::length(accelerations_[0]) + 1e-6f));

    return std::min({dt_cfl, dt_force, params_.timestep});
}

void SPHEngine::update_neighbor_lists() {
    // Build spatial hash
    std::vector<glm::vec3> positions = particles_.get_positions();
    spatial_hash_->build(positions);

    // Query neighbors for each particle
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < particles_.size(); ++i) {
        neighbor_lists_[i] = spatial_hash_->query_squared(
            particles_[i].position,
            params_.neighbor_search_radius * params_.neighbor_search_radius
        );

        perf_stats_.total_neighbor_queries++;
        perf_stats_.max_neighbors = std::max(perf_stats_.max_neighbors, neighbor_lists_[i].size());
    }
}

void SPHEngine::validate_simulation() const {
    if (!utils::validate_particle_setup(std::vector<Particle>(particles_.begin(), particles_.end()))) {
        std::cerr << "Warning: Particle setup validation failed\n";
    }

    float avg_neighbors = utils::compute_average_neighbors(*this);
    if (avg_neighbors < 20.0f || avg_neighbors > 80.0f) {
        std::cerr << "Warning: Average neighbor count (" << avg_neighbors
                  << ") may affect simulation stability\n";
    }
}

// SPH equations implementation
namespace equations {

float compute_density(size_t i, const std::vector<size_t>& neighbors,
                     const ParticleSystem& particles, const Kernel& kernel) {
    const auto& particle_i = particles[i];
    float density = particle_i.mass * kernel.W(glm::vec3(0.0f));  // Self-contribution

    for (size_t j : neighbors) {
        if (i == j) continue;  // Skip self
        const auto& particle_j = particles[j];
        glm::vec3 r_ij = particle_i.position - particle_j.position;
        density += particle_j.mass * kernel.W(r_ij);
    }

    return density;
}

float compute_pressure(float density, float rest_density, float gas_constant) {
    // Ideal gas equation of state
    return gas_constant * (density - rest_density);
}

glm::vec3 compute_pressure_force(size_t i, const std::vector<size_t>& neighbors,
                                const ParticleSystem& particles, const Kernel& kernel,
                                const std::vector<float>& pressures, const std::vector<float>& densities) {
    const auto& particle_i = particles[i];
    glm::vec3 force(0.0f);

    for (size_t j : neighbors) {
        if (i == j) continue;
        const auto& particle_j = particles[j];

        glm::vec3 r_ij = particle_i.position - particle_j.position;
        float r_len = glm::length(r_ij);

        if (r_len < 1e-6f) continue;  // Avoid division by zero

        float pressure_term = (pressures[i] + pressures[j]) / (2.0f * densities[j]);
        glm::vec3 gradW = kernel.gradW(r_ij);

        force -= particle_j.mass * pressure_term * gradW;
    }

    return force;
}

glm::vec3 compute_viscosity_force(size_t i, const std::vector<size_t>& neighbors,
                                 const ParticleSystem& particles, const Kernel& kernel,
                                 float viscosity, const std::vector<float>& densities) {
    const auto& particle_i = particles[i];
    glm::vec3 force(0.0f);

    for (size_t j : neighbors) {
        if (i == j) continue;
        const auto& particle_j = particles[j];

        glm::vec3 r_ij = particle_i.position - particle_j.position;
        glm::vec3 v_ij = particle_j.velocity - particle_i.velocity;

        float laplacianW = kernel.laplacianW(r_ij);
        force += (particle_j.mass / densities[j]) * viscosity * v_ij * laplacianW;
    }

    return force;
}

glm::vec3 compute_external_forces(float gravity) {
    return glm::vec3(0.0f, gravity, 0.0f);
}

glm::vec3 compute_acceleration(const glm::vec3& pressure_force,
                             const glm::vec3& viscosity_force,
                             const glm::vec3& external_force,
                             float particle_mass) {
    return (pressure_force + viscosity_force + external_force * particle_mass) / particle_mass;
}

} // namespace equations

// Utility functions implementation
namespace utils {

std::vector<Particle> create_dam_break_setup(const glm::vec3& dam_size,
                                           const glm::vec3& fluid_size,
                                           float spacing,
                                           const SPHParameters& params) {
    std::vector<Particle> particles;

    // Create dam boundaries first (stationary particles)
    std::vector<Particle> boundaries = create_boundary_box(
        glm::vec3(0.0f, dam_size.y/2.0f, 0.0f),
        dam_size,
        spacing,
        params.particle_mass
    );

    for (auto& particle : boundaries) {
        particle.type = ParticleType::BOUNDARY;
        particle.velocity = glm::vec3(0.0f);  // Stationary
        particle.color = glm::vec3(0.5f, 0.5f, 0.5f);
    }

    // Create fluid particles
    std::vector<Particle> fluid = create_fluid_block(
        glm::vec3(-dam_size.x/2.0f + fluid_size.x/2.0f, fluid_size.y/2.0f, 0.0f),
        fluid_size,
        spacing,
        params.particle_mass
    );

    for (auto& particle : fluid) {
        particle.type = ParticleType::FLUID;
        particle.color = glm::vec3(0.0f, 0.5f, 1.0f);
    }

    particles.insert(particles.end(), boundaries.begin(), boundaries.end());
    particles.insert(particles.end(), fluid.begin(), fluid.end());

    return particles;
}

std::vector<Particle> create_fluid_drop_setup(const glm::vec3& center,
                                            float radius,
                                            float spacing,
                                            const SPHParameters& params) {
    std::vector<Particle> particles;

    int num_particles_per_dim = static_cast<int>(2.0f * radius / spacing);
    glm::vec3 start = center - glm::vec3(radius);

    for (int i = 0; i < num_particles_per_dim; ++i) {
        for (int j = 0; j < num_particles_per_dim; ++j) {
            for (int k = 0; k < num_particles_per_dim; ++k) {
                glm::vec3 position = start + glm::vec3(i * spacing, j * spacing, k * spacing);
                glm::vec3 to_center = position - center;

                if (glm::dot(to_center, to_center) <= radius * radius) {
                    Particle particle(position, params.particle_mass, ParticleType::FLUID);
                    particle.color = glm::vec3(0.0f, 0.7f, 1.0f);
                    particles.push_back(particle);
                }
            }
        }
    }

    return particles;
}

std::vector<Particle> create_granular_flow_setup(const glm::vec3& pile_size,
                                               const glm::vec3& domain_size,
                                               float spacing,
                                               const SPHParameters& params) {
    std::vector<Particle> particles;

    // Create ground and walls
    std::vector<Particle> boundaries = create_boundary_box(
        glm::vec3(0.0f, domain_size.y/2.0f, 0.0f),
        domain_size,
        spacing,
        params.particle_mass
    );

    for (auto& particle : boundaries) {
        particle.type = ParticleType::BOUNDARY;
        particle.color = glm::vec3(0.4f, 0.4f, 0.4f);
    }

    // Create granular pile
    std::vector<Particle> granular = create_fluid_block(
        glm::vec3(0.0f, pile_size.y/2.0f + spacing, 0.0f),
        pile_size,
        spacing,
        params.particle_mass
    );

    for (auto& particle : granular) {
        particle.type = ParticleType::SOLID;
        particle.color = glm::vec3(0.8f, 0.6f, 0.2f);  // Sand color
    }

    particles.insert(particles.end(), boundaries.begin(), boundaries.end());
    particles.insert(particles.end(), granular.begin(), granular.end());

    return particles;
}

bool validate_particle_setup(const std::vector<Particle>& particles) {
    if (particles.empty()) return false;

    // Check for overlapping particles
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            glm::vec3 diff = particles[i].position - particles[j].position;
            if (glm::dot(diff, diff) < 1e-8f) {
                std::cerr << "Warning: Overlapping particles detected\n";
                return false;
            }
        }
    }

    return true;
}

float compute_average_neighbors(const SPHEngine& engine) {
    const auto& particles = engine.get_particles();
    if (particles.size() == 0) return 0.0f;

    // This is a simplified version - in practice you'd need access to neighbor lists
    // For now, return an estimate based on typical SPH neighbor counts
    return 50.0f;  // Typical value for cubic spline kernel
}

} // namespace utils

} // namespace sph
