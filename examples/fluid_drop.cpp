#include "sph_engine.h"
#include "renderer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

/**
 * Fluid Drop Simulation Example
 *
 * This example simulates a spherical drop of fluid falling under gravity,
 * demonstrating surface tension effects and fluid deformation.
 */

int main(int argc, char* argv[]) {
    std::cout << "SPH Fluid Drop Simulation" << std::endl;
    std::cout << "=========================" << std::endl;

    try {
        // Configuration
        int num_particles = 8000;
        bool render = true;
        float drop_radius = 0.15f;

        if (argc > 1) num_particles = std::atoi(argv[1]);
        if (argc > 2) render = std::atoi(argv[2]) != 0;
        if (argc > 3) drop_radius = std::atof(argv[3]);

        // Initialize SPH engine
        sph::SPHEngine engine(num_particles);

        sph::SPHParameters params;
        params.rest_density = 1000.0f;
        params.gas_constant = 2000.0f;
        params.viscosity = 0.01f;  // Higher viscosity for drop stability
        params.smoothing_length = 0.02f;
        params.particle_mass = 0.001f;
        params.timestep = 0.0008f;  // Smaller timestep for stability
        params.gravity = -9.81f;

        // Set simulation domain
        params.bounds.xmin = -0.5f;
        params.bounds.xmax = 0.5f;
        params.bounds.ymin = -0.2f;
        params.bounds.ymax = 0.8f;
        params.bounds.zmin = -0.5f;
        params.bounds.zmax = 0.5f;

        engine.initialize(params);

        // Create fluid drop
        std::vector<sph::Particle> drop_particles = sph::utils::create_fluid_drop_setup(
            glm::vec3(0.0f, 0.4f, 0.0f),  // Position above ground
            drop_radius,
            0.012f,  // Particle spacing
            params
        );

        // Create ground boundary
        std::vector<sph::Particle> ground_particles = sph::create_boundary_box(
            glm::vec3(0.0f, -0.1f, 0.0f),
            glm::vec3(1.0f, 0.2f, 1.0f),
            0.015f,
            params.particle_mass
        );

        engine.add_particles(drop_particles);
        engine.add_particles(ground_particles);

        std::cout << "Simulation initialized:" << std::endl;
        std::cout << "  Fluid particles: " << drop_particles.size() << std::endl;
        std::cout << "  Boundary particles: " << ground_particles.size() << std::endl;
        std::cout << "  Total particles: " << engine.get_particles().size() << std::endl;

        // Initialize renderer
        std::unique_ptr<sph::Renderer> renderer;
        if (render) {
            renderer = std::make_unique<sph::Renderer>();

            sph::RenderParameters render_params;
            render_params.window_width = 1024;
            render_params.window_height = 768;
            render_params.window_title = "SPH Fluid Drop Simulation";
            render_params.particle_size = 0.006f;
            render_params.use_instancing = true;
            render_params.max_instances = num_particles;
            render_params.color_mode = sph::RenderParameters::ColorMode::VELOCITY;
            render_params.background_color = glm::vec3(0.1f, 0.15f, 0.2f);

            if (!renderer->initialize(render_params)) {
                std::cerr << "Failed to initialize renderer" << std::endl;
                render = false;
            } else if (!renderer->load_shaders("shaders/particle.vert", "shaders/particle.frag")) {
                std::cerr << "Failed to load shaders" << std::endl;
                render = false;
            } else {
                // Set camera for drop view
                renderer->get_camera().position = glm::vec3(0.3f, 0.3f, 0.8f);
                renderer->get_camera().yaw = -30.0f;
                renderer->get_camera().pitch = -20.0f;
            }
        }

        // Simulation parameters
        const float simulation_duration = 4.0f;
        const float target_dt = params.timestep;
        const size_t max_steps = static_cast<size_t>(simulation_duration / target_dt);

        // Tracking variables
        float max_spread = 0.0f;
        float min_height = drop_radius;
        float max_height = 0.4f;

        std::cout << std::endl << "Starting simulation..." << std::endl;
        std::cout << std::fixed << std::setprecision(4);

        auto start_time = std::chrono::high_resolution_clock::now();

        for (size_t step = 0; step < max_steps; ++step) {
            // Run simulation step
            engine.step(target_dt);

            // Update renderer
            if (render && renderer) {
                renderer->render_frame(engine);
                if (renderer->should_close()) {
                    break;
                }
            }

            // Analyze drop properties every 50 steps
            if (step % 50 == 0) {
                const auto& positions = engine.get_positions();
                const auto& velocities = engine.get_velocities();

                // Find fluid particles (exclude boundaries)
                float current_max_spread = 0.0f;
                float current_min_height = 1.0f;
                float current_max_height = -1.0f;
                float avg_velocity = 0.0f;

                for (size_t i = 0; i < positions.size(); ++i) {
                    const auto& pos = positions[i];
                    const auto& vel = velocities[i];

                    if (i < drop_particles.size()) {  // Fluid particles
                        float horizontal_dist = std::sqrt(pos.x * pos.x + pos.z * pos.z);
                        current_max_spread = std::max(current_max_spread, horizontal_dist);
                        current_min_height = std::min(current_min_height, pos.y);
                        current_max_height = std::max(current_max_height, pos.y);
                        avg_velocity += glm::length(vel);
                    }
                }

                avg_velocity /= drop_particles.size();
                max_spread = std::max(max_spread, current_max_spread);
                min_height = std::min(min_height, current_min_height);

                // Output progress
                float sim_time = engine.get_current_time();
                std::cout << "t=" << sim_time << "s | Spread: " << current_max_spread
                         << " | Height: " << current_min_height << "-"
                         << current_max_height << " | Avg vel: " << avg_velocity << std::endl;
            }
        }

        // Final analysis
        auto end_time = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(end_time - start_time).count();

        std::cout << std::endl << "Simulation completed!" << std::endl;
        std::cout << "Results:" << std::endl;
        std::cout << "  Maximum spread: " << max_spread << std::endl;
        std::cout << "  Minimum height reached: " << min_height << std::endl;
        std::cout << "  Drop deformation ratio: " << (max_spread / drop_radius) << std::endl;
        std::cout << "  Performance: " << (max_steps / total_time) << " FPS" << std::endl;

        // Conservation check
        float mass_error, energy_error;
        engine.compute_conservation_errors(mass_error, energy_error);
        std::cout << "  Mass conservation error: " << (mass_error * 100.0f) << "%" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

