#include "sph_engine.h"
#include "renderer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

int main(int argc, char* argv[]) {
    std::cout << "SPH Particle Simulator" << std::endl;
    std::cout << "=====================" << std::endl;

    try {
        // Initialize SPH engine
        std::cout << "Initializing SPH Engine..." << std::endl;
        sph::SPHEngine engine(100000);  // Support up to 100k particles

        sph::SPHParameters sph_params;
        sph_params.rest_density = 1000.0f;
        sph_params.gas_constant = 2000.0f;
        sph_params.viscosity = 0.001f;
        sph_params.smoothing_length = 0.02f;
        sph_params.particle_mass = 0.001f;
        sph_params.timestep = 0.001f;
        sph_params.gravity = -9.81f;

        // Set boundaries
        sph_params.bounds.xmin = -1.0f;
        sph_params.bounds.xmax = 1.0f;
        sph_params.bounds.ymin = -1.0f;
        sph_params.bounds.ymax = 1.0f;
        sph_params.bounds.zmin = -1.0f;
        sph_params.bounds.zmax = 1.0f;

        engine.initialize(sph_params);

        // Initialize dam break scenario
        std::cout << "Setting up dam break simulation..." << std::endl;
        engine.initialize_dam_break();

        std::cout << "Initialized with " << engine.get_particles().size() << " particles" << std::endl;

        // Initialize renderer
        std::cout << "Initializing OpenGL Renderer..." << std::endl;
        sph::Renderer renderer;

        sph::RenderParameters render_params;
        render_params.window_width = 1280;
        render_params.window_height = 720;
        render_params.window_title = "SPH Particle Simulator - Dam Break";
        render_params.particle_size = 0.005f;
        render_params.use_instancing = true;
        render_params.max_instances = 100000;
        render_params.color_mode = sph::RenderParameters::ColorMode::VELOCITY;
        render_params.background_color = glm::vec3(0.1f, 0.1f, 0.15f);

        if (!renderer.initialize(render_params)) {
            std::cerr << "Failed to initialize renderer" << std::endl;
            return -1;
        }

        // Load shaders
        std::cout << "Loading shaders..." << std::endl;
        if (!renderer.load_shaders("shaders/particle.vert", "shaders/particle.frag")) {
            std::cerr << "Failed to load shaders" << std::endl;
            return -1;
        }

        // Set camera position for dam break view
        renderer.get_camera().position = glm::vec3(0.0f, 0.0f, 2.0f);
        renderer.get_camera().yaw = -90.0f;
        renderer.get_camera().pitch = 0.0f;

        // Simulation loop
        std::cout << "Starting simulation loop..." << std::endl;
        std::cout << "Controls:" << std::endl;
        std::cout << "  WASD: Move camera" << std::endl;
        std::cout << "  Mouse: Look around (right-click to capture)" << std::endl;
        std::cout << "  ESC: Release mouse" << std::endl;
        std::cout << "  Q/E: Move up/down" << std::endl;
        std::cout << std::endl;

        size_t frame_count = 0;
        auto start_time = std::chrono::high_resolution_clock::now();

        while (!renderer.should_close()) {
            // Run simulation step
            engine.step();  // Adaptive timestep

            // Render frame
            renderer.render_frame(engine);

            // Print stats every 100 frames
            frame_count++;
            if (frame_count % 100 == 0) {
                auto current_time = std::chrono::high_resolution_clock::now();
                double elapsed = std::chrono::duration<double>(current_time - start_time).count();
                double fps = frame_count / elapsed;

                std::cout << std::fixed << std::setprecision(2);
                std::cout << "Frame: " << frame_count
                         << " | Sim Time: " << engine.get_current_time()
                         << " | FPS: " << fps
                         << " | Particles: " << engine.get_particles().size()
                         << " | Render FPS: " << renderer.get_fps()
                         << std::endl;

                // Conservation check
                float mass_error, energy_error;
                engine.compute_conservation_errors(mass_error, energy_error);
                std::cout << "Mass conservation error: " << (mass_error * 100.0f) << "%" << std::endl;
            }

            // Limit frame rate if needed
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

        std::cout << "Simulation finished. Total frames: " << frame_count << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

