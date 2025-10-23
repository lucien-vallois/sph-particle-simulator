#include "renderer.h"
#include "sph_engine.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <glm/gtc/type_ptr.hpp>

namespace sph {

// Camera implementation
void Camera::update_vectors() {
    glm::vec3 front_vec;
    front_vec.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front_vec.y = sin(glm::radians(pitch));
    front_vec.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    front = glm::normalize(front_vec);

    right = glm::normalize(glm::cross(front, glm::vec3(0.0f, 1.0f, 0.0f)));
    up = glm::normalize(glm::cross(right, front));
}

glm::mat4 Camera::get_view_matrix() const {
    return glm::lookAt(position, position + front, up);
}

glm::mat4 Camera::get_projection_matrix(float aspect_ratio) const {
    return glm::perspective(glm::radians(fov), aspect_ratio, 0.1f, 100.0f);
}

// Renderer implementation
Renderer::Renderer() = default;

Renderer::~Renderer() {
    if (VAO_) glDeleteVertexArrays(1, &VAO_);
    if (VBO_) glDeleteBuffers(1, &VBO_);
    if (instance_VBO_) glDeleteBuffers(1, &instance_VBO_);
    if (shader_program_) glDeleteProgram(shader_program_);
    if (window_) glfwDestroyWindow(window_);

    glfwTerminate();
}

bool Renderer::initialize(const RenderParameters& params) {
    params_ = params;

    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return false;
    }

    // Configure GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);  // MSAA

    // Create window
    window_ = glfwCreateWindow(params_.window_width, params_.window_height,
                              params_.window_title.c_str(), nullptr, nullptr);
    if (!window_) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(window_);
    glfwSetWindowUserPointer(window_, this);

    // Set callbacks
    glfwSetFramebufferSizeCallback(window_, framebuffer_size_callback);
    glfwSetCursorPosCallback(window_, mouse_callback);
    glfwSetScrollCallback(window_, scroll_callback);

    // Configure OpenGL
    if (!setup_opengl()) {
        return false;
    }

    // Create vertex buffers
    if (!create_vertex_buffers()) {
        return false;
    }

    // Setup instancing
    setup_instancing(params_.max_instances);

    // Enable depth testing and blending
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Set background color
    glClearColor(params_.background_color.r, params_.background_color.g,
                 params_.background_color.b, 1.0f);

    return true;
}

bool Renderer::setup_opengl() {
    // Load OpenGL functions (assuming GLAD or similar is available)
    // In a real implementation, you'd initialize GLAD here
    // For now, we'll assume the functions are available

    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    return true;
}

bool Renderer::create_vertex_buffers() {
    // Create VAO
    glGenVertexArrays(1, &VAO_);
    glBindVertexArray(VAO_);

    // Create VBO for quad vertices (for instanced rendering)
    float quad_vertices[] = {
        // positions
        -0.5f, -0.5f, 0.0f,
         0.5f, -0.5f, 0.0f,
         0.5f,  0.5f, 0.0f,
        -0.5f,  0.5f, 0.0f,
    };

    glGenBuffers(1, &VBO_);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad_vertices), quad_vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Create instance VBO for particle data
    glGenBuffers(1, &instance_VBO_);
    glBindBuffer(GL_ARRAY_BUFFER, instance_VBO_);
    // Buffer will be filled dynamically

    // Position (location 1)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribDivisor(1, 1);

    // Velocity (location 2) - for color mapping
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1);

    // Color (location 3)
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1);

    glBindVertexArray(0);
    return true;
}

void Renderer::setup_instancing(size_t max_particles) {
    // Allocate buffer for instance data
    glBindBuffer(GL_ARRAY_BUFFER, instance_VBO_);
    glBufferData(GL_ARRAY_BUFFER, max_particles * 9 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
}

bool Renderer::load_shaders(const std::string& vertex_path, const std::string& fragment_path) {
    std::string vertex_source = load_shader_source(vertex_path);
    std::string fragment_source = load_shader_source(fragment_path);

    if (vertex_source.empty() || fragment_source.empty()) {
        return false;
    }

    unsigned int vertex_shader = compile_shader(vertex_source, GL_VERTEX_SHADER);
    unsigned int fragment_shader = compile_shader(fragment_source, GL_FRAGMENT_SHADER);

    if (vertex_shader == 0 || fragment_shader == 0) {
        return false;
    }

    if (!link_program(vertex_shader, fragment_shader)) {
        return false;
    }

    // Get uniform locations
    view_loc_ = glGetUniformLocation(shader_program_, "view");
    projection_loc_ = glGetUniformLocation(shader_program_, "projection");
    model_loc_ = glGetUniformLocation(shader_program_, "model");
    particle_size_loc_ = glGetUniformLocation(shader_program_, "particleSize");
    color_mode_loc_ = glGetUniformLocation(shader_program_, "colorMode");

    return true;
}

unsigned int Renderer::compile_shader(const std::string& source, unsigned int type) {
    unsigned int shader = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    int success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char info_log[512];
        glGetShaderInfoLog(shader, 512, nullptr, info_log);
        std::cerr << "Shader compilation failed:\n" << info_log << std::endl;
        glDeleteShader(shader);
        return 0;
    }

    return shader;
}

bool Renderer::link_program(unsigned int vertex_shader, unsigned int fragment_shader) {
    shader_program_ = glCreateProgram();
    glAttachShader(shader_program_, vertex_shader);
    glAttachShader(shader_program_, fragment_shader);
    glLinkProgram(shader_program_);

    int success;
    glGetProgramiv(shader_program_, GL_LINK_STATUS, &success);
    if (!success) {
        char info_log[512];
        glGetProgramInfoLog(shader_program_, 512, nullptr, info_log);
        std::cerr << "Program linking failed:\n" << info_log << std::endl;
        return false;
    }

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    return true;
}

std::string Renderer::load_shader_source(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open shader file: " << filepath << std::endl;
        return "";
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

void Renderer::render_frame(const SPHEngine& engine) {
    // Update timing
    double current_time = glfwGetTime();
    delta_time_ = current_time - last_frame_time_;
    last_frame_time_ = current_time;

    update_fps_counter();

    // Process input
    process_input();

    // Clear framebuffer
    clear_framebuffer();

    // Update camera
    camera_.update_vectors();

    // Update particle data
    update_particle_data(engine);

    // Render particles
    if (params_.use_instancing) {
        render_particles_instanced(engine.get_particles().size());
    } else {
        render_particles_individual(engine.get_particles().size());
    }

    // Swap buffers
    glfwSwapBuffers(window_);
    glfwPollEvents();
}

void Renderer::clear_framebuffer() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void Renderer::update_particle_data(const SPHEngine& engine) {
    const auto& particles = engine.get_particles();
    const auto& positions = engine.get_positions();
    const auto& velocities = engine.get_velocities();

    if (particles.size() == 0) return;

    // Prepare instance data: position (3), velocity (3), color (3)
    std::vector<float> instance_data;
    instance_data.reserve(particles.size() * 9);

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];

        // Position
        instance_data.push_back(particle.position.x);
        instance_data.push_back(particle.position.y);
        instance_data.push_back(particle.position.z);

        // Velocity (for color mapping)
        instance_data.push_back(particle.velocity.x);
        instance_data.push_back(particle.velocity.y);
        instance_data.push_back(particle.velocity.z);

        // Color
        instance_data.push_back(particle.color.r);
        instance_data.push_back(particle.color.g);
        instance_data.push_back(particle.color.b);
    }

    // Update instance buffer
    glBindBuffer(GL_ARRAY_BUFFER, instance_VBO_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, instance_data.size() * sizeof(float), instance_data.data());
}

void Renderer::render_particles_instanced(size_t particle_count) {
    glUseProgram(shader_program_);
    glBindVertexArray(VAO_);

    // Set uniforms
    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    float aspect_ratio = static_cast<float>(width) / height;

    glUniformMatrix4fv(view_loc_, 1, GL_FALSE, glm::value_ptr(camera_.get_view_matrix()));
    glUniformMatrix4fv(projection_loc_, 1, GL_FALSE, glm::value_ptr(camera_.get_projection_matrix(aspect_ratio)));
    glUniformMatrix4fv(model_loc_, 1, GL_FALSE, glm::value_ptr(glm::mat4(1.0f)));
    glUniform1f(particle_size_loc_, params_.particle_size);
    glUniform1i(color_mode_loc_, static_cast<int>(params_.color_mode));

    // Render instanced quads
    glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, static_cast<GLsizei>(particle_count));

    glBindVertexArray(0);
}

void Renderer::render_particles_individual(size_t particle_count) {
    // Fallback for non-instanced rendering (slower)
    glUseProgram(shader_program_);

    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    float aspect_ratio = static_cast<float>(width) / height;

    glUniformMatrix4fv(view_loc_, 1, GL_FALSE, glm::value_ptr(camera_.get_view_matrix()));
    glUniformMatrix4fv(projection_loc_, 1, GL_FALSE, glm::value_ptr(camera_.get_projection_matrix(aspect_ratio)));
    glUniform1f(particle_size_loc_, params_.particle_size);
    glUniform1i(color_mode_loc_, static_cast<int>(params_.color_mode));

    // Individual rendering (much slower for large particle counts)
    glBindVertexArray(VAO_);
    for (size_t i = 0; i < particle_count; ++i) {
        // Would need individual model matrices - not implemented for performance
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    }
    glBindVertexArray(0);
}

bool Renderer::should_close() const {
    return glfwWindowShouldClose(window_);
}

void Renderer::process_input() {
    process_keyboard_input();

    // Mouse capture toggle
    if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        mouse_captured_ = false;
        glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    if (glfwGetMouseButton(window_, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        mouse_captured_ = true;
        glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        first_mouse_ = true;
    }
}

void Renderer::process_keyboard_input() {
    float camera_speed = 2.5f * static_cast<float>(delta_time_);

    if (glfwGetKey(window_, GLFW_KEY_W) == GLFW_PRESS)
        camera_.position += camera_speed * camera_.front;
    if (glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS)
        camera_.position -= camera_speed * camera_.front;
    if (glfwGetKey(window_, GLFW_KEY_A) == GLFW_PRESS)
        camera_.position -= camera_speed * camera_.right;
    if (glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS)
        camera_.position += camera_speed * camera_.right;
    if (glfwGetKey(window_, GLFW_KEY_Q) == GLFW_PRESS)
        camera_.position -= camera_speed * camera_.up;
    if (glfwGetKey(window_, GLFW_KEY_E) == GLFW_PRESS)
        camera_.position += camera_speed * camera_.up;
}

// Static callback functions
void Renderer::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

void Renderer::mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));

    if (!renderer->mouse_captured_) return;

    if (renderer->first_mouse_) {
        renderer->last_mouse_x_ = xpos;
        renderer->last_mouse_y_ = ypos;
        renderer->first_mouse_ = false;
    }

    float xoffset = static_cast<float>(xpos - renderer->last_mouse_x_);
    float yoffset = static_cast<float>(renderer->last_mouse_y_ - ypos);

    renderer->last_mouse_x_ = xpos;
    renderer->last_mouse_y_ = ypos;

    const float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    renderer->camera_.yaw += xoffset;
    renderer->camera_.pitch += yoffset;

    renderer->camera_.pitch = std::clamp(renderer->camera_.pitch, -89.0f, 89.0f);
}

void Renderer::scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
    renderer->camera_.fov -= static_cast<float>(yoffset);
    renderer->camera_.fov = std::clamp(renderer->camera_.fov, 1.0f, 45.0f);
}

void Renderer::update_fps_counter() {
    frame_count_++;
    if (current_time_ - fps_update_time_ >= 1.0) {
        current_fps_ = frame_count_ / (current_time_ - fps_update_time_);
        frame_count_ = 0;
        fps_update_time_ = current_time_;

        // Update window title with FPS
        std::string title = params_.window_title + " - FPS: " + std::to_string(static_cast<int>(current_fps_));
        glfwSetWindowTitle(window_, title.c_str());
    }
}

void Renderer::set_particle_size(float size) {
    params_.particle_size = size;
}

void Renderer::set_color_mode(RenderParameters::ColorMode mode) {
    params_.color_mode = mode;
}

void Renderer::set_background_color(const glm::vec3& color) {
    params_.background_color = color;
    glClearColor(color.r, color.g, color.b, 1.0f);
}

void Renderer::get_window_size(int& width, int& height) const {
    glfwGetWindowSize(window_, &width, &height);
}

// Utility functions
glm::vec3 Renderer::hsv_to_rgb(float h, float s, float v) {
    float c = v * s;
    float x = c * (1.0f - std::abs(std::fmod(h / 60.0f, 2.0f) - 1.0f));
    float m = v - c;

    float r, g, b;
    if (h >= 0 && h < 60) {
        r = c; g = x; b = 0;
    } else if (h >= 60 && h < 120) {
        r = x; g = c; b = 0;
    } else if (h >= 120 && h < 180) {
        r = 0; g = c; b = x;
    } else if (h >= 180 && h < 240) {
        r = 0; g = x; b = c;
    } else if (h >= 240 && h < 300) {
        r = x; g = 0; b = c;
    } else {
        r = c; g = 0; b = x;
    }

    return glm::vec3(r + m, g + m, b + m);
}

float Renderer::map_to_color(float value, float min_val, float max_val) {
    if (max_val == min_val) return 0.5f;
    return (value - min_val) / (max_val - min_val);
}

} // namespace sph

