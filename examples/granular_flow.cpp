#include "sph_engine.h"
#include "renderer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

/**
 * Granular Flow Simulation Example
 *
 * This example simulates granular material flow, useful for modeling
 * sand, soil, or other granular materials in civil engineering and
 * geotechnical applications.
 */

int main(int argc, char* argv[]) {
    std::cout << "SPH Granular Flow Simulation" << std::endl;
    std::cout << "============================" << std::endl;

    try {
        // Configuration
        int num_particles = 12000;
        bool render = true;

        if (argc > 1) num_particles = std::atoi(argv[1]);
        if (argc > 2) render = std::atoi(argv[2]) != 0;

        // Initialize SPH engine
        sph::SPHEngine engine(num_particles);

        sph::SPHParameters params;
        params.rest_density = 1500.0f;     // Higher density for granular material
        params.gas_constant = 5000.0f;     // Higher pressure for granular behavior
        params.viscosity = 0.05f;          // Higher viscosity for granular friction
        params.smoothing_length = 0.022f;
        params.particle_mass = 0.002f;     // Heavier particles
        params.timestep = 0.0005f;         // Smaller timestep for stability
        params.gravity = -9.81f;
        params.damping = 0.98f;            // Higher damping for energy dissipation

        // Set simulation domain
        params.bounds.xmin = -0.8f;
        params.bounds.xmax = 0.8f;
        params.bounds.ymin = -0.3f;
        params.bounds.ymax = 0.7f;
        params.bounds.zmin = -0.8f;
        params.bounds.zmax = 0.8f;

        engine.initialize(params);
        engine.initialize_granular_flow();

        std::cout << "Simulation initialized with " << engine.get_particles().size() << " particles" << std::endl;

        // Initialize renderer
        std::unique_ptr<sph::Renderer> renderer;
        if (render) {
            renderer = std::make_unique<sph::Renderer>();

            sph::RenderParameters render_params;
            render_params.window_width = 1200;
            render_params.window_height = 800;
            render_params.window_title = "SPH Granular Flow Simulation";
            render_params.particle_size = 0.008f;
            render_params.use_instancing = true;
            render_params.max_instances = num_particles;
            render_params.color_mode = sph::RenderParameters::ColorMode::TYPE;  // Color by particle type
            render_params.background_color = glm::vec3(0.15f, 0.12f, 0.1f);

            if (!renderer->initialize(render_params)) {
                std::cerr << "Failed to initialize renderer" << std::endl;
                render = false;
            } else if (!renderer->load_shaders("shaders/particle.vert", "shaders/particle.frag")) {
                std::cerr << "Failed to load shaders" << std::endl;
                render = false;
            } else {
                // Set camera for flow view
                renderer->get_camera().position = glm::vec3(1.0f, 0.8f, 1.5f);
                renderer->get_camera().yaw = -35.0f;
                renderer->get_camera().pitch = -25.0f;
            }
        }

        // Simulation parameters
        const float simulation_duration = 5.0f;
        const float target_dt = params.timestep;
        const size_t max_steps = static_cast<size_t>(simulation_duration / target_dt);

        // Flow analysis variables
        float initial_pile_height = 0.4f;
        float flow_distance = 0.0f;
        float max_flow_velocity = 0.0f;
        std::vector<float> flow_front_positions;

        std::cout << std::endl << "Starting granular flow simulation..." << std::endl;
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

            // Analyze flow properties every 100 steps
            if (step % 100 == 0) {
                const auto& positions = engine.get_positions();
                const auto& velocities = engine.get_velocities();

                // Find granular particles (exclude boundaries)
                float min_x = 1.0f;
                float max_x = -1.0f;
                float avg_velocity = 0.0f;
                int moving_particles = 0;

                for (size_t i = 0; i < positions.size(); ++i) {
                    const auto& pos = positions[i];
                    const auto& vel = velocities[i];

                    // Assume first particles are granular (from initialize_granular_flow)
                    if (i < num_particles * 0.8f) {  // Roughly granular particles
                        min_x = std::min(min_x, pos.x);
                        max_x = std::max(max_x, pos.x);

                        float speed = glm::length(vel);
                        avg_velocity += speed;
                        max_flow_velocity = std::max(max_flow_velocity, speed);

                        if (speed > 0.01f) {  // Threshold for "moving"
                            moving_particles++;
                        }
                    }
                }

                avg_velocity /= (num_particles * 0.8f);
                flow_distance = std::max(flow_distance, std::abs(max_x) + std::abs(min_x));

                // Track flow front
                flow_front_positions.push_back(max_x);

                // Output progress
                float sim_time = engine.get_current_time();
                std::cout << "t=" << sim_time << "s | Flow extent: " << min_x
                         << " to " << max_x << " | Distance: " << flow_distance
                         << " | Moving particles: " << moving_particles
                         << " | Avg vel: " << avg_velocity << std::endl;
            }
        }

        // Final analysis
        auto end_time = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(end_time - start_time).count();

        std::cout << std::endl << "Granular flow simulation completed!" << std::endl;
        std::cout << "Flow characteristics:" << std::endl;
        std::cout << "  Maximum flow distance: " << flow_distance << std::endl;
        std::cout << "  Maximum particle velocity: " << max_flow_velocity << std::endl;
        std::cout << "  Final flow front position: " << flow_front_positions.back() << std::endl;

        // Calculate flow mobility
        float initial_width = 0.3f;  // From initialize_granular_flow
        float flow_mobility = flow_distance / initial_width;
        std::cout << "  Flow mobility (distance/width): " << flow_mobility << std::endl;

        // Performance
        std::cout << "Performance: " << (max_steps / total_time) << " FPS" << std::endl;

        // Conservation
        float mass_error, energy_error;
        engine.compute_conservation_errors(mass_error, energy_error);
        std::cout << "Mass conservation error: " << (mass_error * 100.0f) << "%" << std::endl;

        // Granular flow analysis
        if (!flow_front_positions.empty()) {
            // Calculate flow velocity (simplified)
            float flow_velocity = 0.0f;
            if (flow_front_positions.size() > 1) {
                float dt = simulation_duration / flow_front_positions.size();
                flow_velocity = (flow_front_positions.back() - flow_front_positions.front()) / simulation_duration;
            }
            std::cout << "Average flow velocity: " << flow_velocity << " m/s" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

