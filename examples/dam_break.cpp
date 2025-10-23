#include "sph_engine.h"
#include "renderer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <fstream>

/**
 * Dam Break Simulation Example
 *
 * This example demonstrates a classic SPH dam break scenario where a column
 * of fluid is released and flows due to gravity, creating complex fluid dynamics.
 */

int main(int argc, char* argv[]) {
    std::cout << "SPH Dam Break Simulation" << std::endl;
    std::cout << "========================" << std::endl;

    try {
        // Parse command line arguments
        int num_particles = 10000;
        bool render = true;
        bool save_data = false;
        float duration = 3.0f;

        if (argc > 1) num_particles = std::atoi(argv[1]);
        if (argc > 2) render = std::atoi(argv[2]) != 0;
        if (argc > 3) save_data = std::atoi(argv[3]) != 0;
        if (argc > 4) duration = std::atof(argv[4]);

        std::cout << "Configuration:" << std::endl;
        std::cout << "  Particles: " << num_particles << std::endl;
        std::cout << "  Rendering: " << (render ? "Enabled" : "Disabled") << std::endl;
        std::cout << "  Data export: " << (save_data ? "Enabled" : "Disabled") << std::endl;
        std::cout << "  Duration: " << duration << " seconds" << std::endl;
        std::cout << std::endl;

        // Initialize SPH engine
        sph::SPHEngine engine(num_particles * 2);  // Extra space for boundaries

        sph::SPHParameters params;
        params.rest_density = 1000.0f;
        params.gas_constant = 2000.0f;
        params.viscosity = 0.001f;
        params.smoothing_length = 0.025f;
        params.particle_mass = 0.001f;
        params.timestep = 0.001f;
        params.gravity = -9.81f;
        params.damping = 0.995f;  // Slight damping for stability

        // Set simulation domain
        params.bounds.xmin = -1.0f;
        params.bounds.xmax = 1.0f;
        params.bounds.ymin = -0.5f;
        params.bounds.ymax = 1.5f;
        params.bounds.zmin = -1.0f;
        params.bounds.zmax = 1.0f;

        engine.initialize(params);
        engine.initialize_dam_break();

        std::cout << "Simulation initialized with " << engine.get_particles().size() << " particles" << std::endl;

        // Initialize renderer if requested
        std::unique_ptr<sph::Renderer> renderer;
        if (render) {
            renderer = std::make_unique<sph::Renderer>();

            sph::RenderParameters render_params;
            render_params.window_width = 1280;
            render_params.window_height = 720;
            render_params.window_title = "SPH Dam Break Simulation";
            render_params.particle_size = 0.008f;
            render_params.use_instancing = true;
            render_params.max_instances = num_particles * 2;
            render_params.color_mode = sph::RenderParameters::ColorMode::VELOCITY;
            render_params.background_color = glm::vec3(0.05f, 0.05f, 0.1f);

            if (!renderer->initialize(render_params)) {
                std::cerr << "Failed to initialize renderer, continuing without visualization" << std::endl;
                renderer.reset();
                render = false;
            } else if (!renderer->load_shaders("shaders/particle.vert", "shaders/particle.frag")) {
                std::cerr << "Failed to load shaders, continuing without visualization" << std::endl;
                renderer.reset();
                render = false;
            } else {
                // Set camera for good view of dam break
                renderer->get_camera().position = glm::vec3(0.5f, 0.5f, 2.0f);
                renderer->get_camera().yaw = -45.0f;
                renderer->get_camera().pitch = -15.0f;
            }
        }

        // Data export setup
        std::ofstream data_file;
        if (save_data) {
            data_file.open("dam_break_data.csv");
            data_file << "time,particles,mass_error,energy,total_energy,avg_density,max_velocity\n";
        }

        // Simulation loop
        const float target_dt = params.timestep;
        const size_t max_steps = static_cast<size_t>(duration / target_dt);
        size_t step_count = 0;

        auto start_time = std::chrono::high_resolution_clock::now();
        auto last_output_time = start_time;

        std::cout << "Starting simulation..." << std::endl;
        std::cout << std::fixed << std::setprecision(3);

        while (step_count < max_steps) {
            // Run simulation step
            engine.step(target_dt);

            // Update renderer
            if (render && renderer) {
                renderer->render_frame(engine);
                if (renderer->should_close()) {
                    break;
                }
            }

            step_count++;

            // Periodic output
            auto current_time = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(current_time - last_output_time).count();

            if (elapsed >= 0.1) {  // Update every 100ms
                double sim_time = engine.get_current_time();
                double real_time_fps = step_count / std::chrono::duration<double>(current_time - start_time).count();

                // Compute statistics
                float mass_error, energy_error;
                engine.compute_conservation_errors(mass_error, energy_error);

                const auto& velocities = engine.get_velocities();
                float max_velocity = 0.0f;
                for (const auto& vel : velocities) {
                    max_velocity = std::max(max_velocity, glm::length(vel));
                }

                const auto& densities = engine.get_densities();
                float avg_density = 0.0f;
                if (!densities.empty()) {
                    for (float d : densities) avg_density += d;
                    avg_density /= densities.size();
                }

                // Console output
                std::cout << "Time: " << sim_time << "s | Step: " << step_count
                         << " | FPS: " << real_time_fps
                         << " | Mass err: " << (mass_error * 100.0f) << "%"
                         << " | Max vel: " << max_velocity
                         << " | Avg density: " << avg_density << std::endl;

                // Data export
                if (save_data && data_file.is_open()) {
                    data_file << sim_time << "," << engine.get_particles().size() << ","
                             << mass_error << "," << energy_error << ","
                             << engine.get_total_energy() << "," << avg_density << ","
                             << max_velocity << "\n";
                }

                last_output_time = current_time;
            }

            // Small delay to prevent excessive CPU usage when not rendering
            if (!render) {
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            }
        }

        // Final statistics
        auto end_time = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(end_time - start_time).count();

        std::cout << std::endl << "Simulation completed!" << std::endl;
        std::cout << "Total simulation time: " << engine.get_current_time() << " seconds" << std::endl;
        std::cout << "Total real time: " << total_time << " seconds" << std::endl;
        std::cout << "Average performance: " << (step_count / total_time) << " FPS" << std::endl;

        // Performance breakdown
        const auto& stats = engine.get_performance_stats();
        if (stats.total_time > 0.0) {
            double total_sim_time = stats.total_time;
            std::cout << std::endl << "Performance breakdown:" << std::endl;
            std::cout << "  Neighbor search: " << (stats.neighbor_search_time / total_sim_time * 100.0) << "%" << std::endl;
            std::cout << "  Density computation: " << (stats.density_computation_time / total_sim_time * 100.0) << "%" << std::endl;
            std::cout << "  Force computation: " << (stats.force_computation_time / total_sim_time * 100.0) << "%" << std::endl;
            std::cout << "  Integration: " << (stats.integration_time / total_sim_time * 100.0) << "%" << std::endl;
        }

        // Conservation properties
        float final_mass_error, final_energy_error;
        engine.compute_conservation_errors(final_mass_error, final_energy_error);
        std::cout << std::endl << "Conservation properties:" << std::endl;
        std::cout << "  Final mass conservation error: " << (final_mass_error * 100.0) << "%" << std::endl;
        std::cout << "  Final energy conservation error: " << (final_energy_error * 100.0) << "%" << std::endl;

        // Cleanup
        if (data_file.is_open()) {
            data_file.close();
            std::cout << "Data exported to 'dam_break_data.csv'" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

