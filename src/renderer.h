#pragma once

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <memory>
#include <string>

namespace sph {

// Forward declarations
class SPHEngine;

// Camera for 3D navigation
struct Camera {
    glm::vec3 position = glm::vec3(0.0f, 0.0f, 3.0f);
    glm::vec3 front = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 right = glm::vec3(1.0f, 0.0f, 0.0f);

    float yaw = -90.0f;
    float pitch = 0.0f;
    float fov = 45.0f;

    void update_vectors();
    glm::mat4 get_view_matrix() const;
    glm::mat4 get_projection_matrix(float aspect_ratio) const;
};

// Rendering parameters
struct RenderParameters {
    // Window settings
    int window_width = 1280;
    int window_height = 720;
    std::string window_title = "SPH Particle Simulator";

    // Visual settings
    float particle_size = 0.01f;
    bool use_instancing = true;
    int max_instances = 1000000;

    // Color mapping
    enum class ColorMode {
        VELOCITY,
        DENSITY,
        PRESSURE,
        TYPE
    } color_mode = ColorMode::VELOCITY;

    // Background
    glm::vec3 background_color = glm::vec3(0.1f, 0.1f, 0.15f);

    // Performance
    bool vsync = true;
    int target_fps = 60;
};

// OpenGL Renderer class
class Renderer {
private:
    // GLFW window
    GLFWwindow* window_ = nullptr;
    RenderParameters params_;

    // OpenGL objects
    unsigned int VAO_ = 0;
    unsigned int VBO_ = 0;
    unsigned int instance_VBO_ = 0;
    unsigned int shader_program_ = 0;

    // Shader locations
    int view_loc_ = -1;
    int projection_loc_ = -1;
    int model_loc_ = -1;
    int particle_size_loc_ = -1;
    int color_mode_loc_ = -1;

    // Camera
    Camera camera_;

    // Performance tracking
    double last_frame_time_ = 0.0;
    double delta_time_ = 0.0;
    int frame_count_ = 0;
    double fps_update_time_ = 0.0;
    double current_fps_ = 0.0;

    // Input handling
    bool first_mouse_ = true;
    double last_mouse_x_ = 0.0;
    double last_mouse_y_ = 0.0;
    bool mouse_captured_ = false;

public:
    // Constructor/Destructor
    Renderer();
    ~Renderer();

    // Initialization
    bool initialize(const RenderParameters& params);
    bool load_shaders(const std::string& vertex_path, const std::string& fragment_path);

    // Main rendering loop
    void render_frame(const SPHEngine& engine);
    bool should_close() const;
    void process_input();

    // Camera control
    void set_camera_position(const glm::vec3& position);
    void set_camera_target(const glm::vec3& target);
    Camera& get_camera() { return camera_; }

    // Rendering parameters
    void set_particle_size(float size);
    void set_color_mode(RenderParameters::ColorMode mode);
    void set_background_color(const glm::vec3& color);

    // Utility
    void get_window_size(int& width, int& height) const;
    double get_delta_time() const { return delta_time_; }
    double get_fps() const { return current_fps_; }

private:
    // OpenGL setup
    bool setup_opengl();
    bool create_vertex_buffers();
    void setup_instancing(size_t max_particles);

    // Shader management
    unsigned int compile_shader(const std::string& source, unsigned int type);
    bool link_program(unsigned int vertex_shader, unsigned int fragment_shader);
    std::string load_shader_source(const std::string& filepath);

    // Rendering functions
    void clear_framebuffer();
    void update_particle_data(const SPHEngine& engine);
    void render_particles_instanced(size_t particle_count);
    void render_particles_individual(size_t particle_count);

    // Input callbacks
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    static void mouse_callback(GLFWwindow* window, double xpos, double ypos);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
    void process_keyboard_input();

    // Utility
    void update_fps_counter();
    glm::vec3 hsv_to_rgb(float h, float s, float v);
    float map_to_color(float value, float min_val, float max_val);
};

// Utility functions for color mapping
namespace render_utils {

// Color mapping functions
glm::vec3 velocity_to_color(const glm::vec3& velocity, float max_velocity = 2.0f);
glm::vec3 density_to_color(float density, float rest_density, float max_density);
glm::vec3 pressure_to_color(float pressure, float max_pressure = 10000.0f);
glm::vec3 particle_type_to_color(ParticleType type);

// Compute value ranges for color mapping
void compute_value_ranges(const SPHEngine& engine,
                         float& min_velocity, float& max_velocity,
                         float& min_density, float& max_density,
                         float& min_pressure, float& max_pressure);

} // namespace render_utils

} // namespace sph

