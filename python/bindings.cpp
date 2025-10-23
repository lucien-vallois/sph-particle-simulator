#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "sph_engine.h"

namespace py = pybind11;

PYBIND11_MODULE(sph, m) {
    m.doc() = "SPH Particle Simulator - Python bindings";

    // SPH Parameters
    py::class_<sph::SPHParameters>(m, "SPHParameters")
        .def(py::init<>())
        .def_readwrite("rest_density", &sph::SPHParameters::rest_density)
        .def_readwrite("gas_constant", &sph::SPHParameters::gas_constant)
        .def_readwrite("viscosity", &sph::SPHParameters::viscosity)
        .def_readwrite("smoothing_length", &sph::SPHParameters::smoothing_length)
        .def_readwrite("particle_mass", &sph::SPHParameters::particle_mass)
        .def_readwrite("timestep", &sph::SPHParameters::timestep)
        .def_readwrite("gravity", &sph::SPHParameters::gravity)
        .def_readwrite("damping", &sph::SPHParameters::damping);

    // Particle type enum
    py::enum_<sph::ParticleType>(m, "ParticleType")
        .value("FLUID", sph::ParticleType::FLUID)
        .value("BOUNDARY", sph::ParticleType::BOUNDARY)
        .value("SOLID", sph::ParticleType::SOLID);

    // Particle class
    py::class_<sph::Particle>(m, "Particle")
        .def(py::init<>())
        .def(py::init<const glm::vec3&, float, sph::ParticleType>(),
             py::arg("position"), py::arg("mass") = 1.0f,
             py::arg("type") = sph::ParticleType::FLUID)
        .def_readwrite("position", &sph::Particle::position)
        .def_readwrite("velocity", &sph::Particle::velocity)
        .def_readwrite("acceleration", &sph::Particle::acceleration)
        .def_readwrite("density", &sph::Particle::density)
        .def_readwrite("pressure", &sph::Particle::pressure)
        .def_readwrite("mass", &sph::Particle::mass)
        .def_readwrite("type", &sph::Particle::type)
        .def_readwrite("temperature", &sph::Particle::temperature)
        .def_readwrite("viscosity", &sph::Particle::viscosity)
        .def_readwrite("color", &sph::Particle::color);

    // Particle System (read-only access to particles)
    py::class_<sph::ParticleSystem>(m, "ParticleSystem")
        .def("size", &sph::ParticleSystem::size)
        .def("capacity", &sph::ParticleSystem::capacity)
        .def("empty", &sph::ParticleSystem::empty)
        .def("get_positions",
             [](const sph::ParticleSystem& ps) {
                 return py::array_t<float>({ps.size(), 3UL},
                     {3 * sizeof(float), sizeof(float)},
                     reinterpret_cast<const float*>(ps.get_positions().data()));
             })
        .def("get_velocities",
             [](const sph::ParticleSystem& ps) {
                 return py::array_t<float>({ps.size(), 3UL},
                     {3 * sizeof(float), sizeof(float)},
                     reinterpret_cast<const float*>(ps.get_velocities().data()));
             })
        .def("get_densities",
             [](const sph::ParticleSystem& ps) {
                 return py::array_t<float>(ps.size(), ps.get_densities().data());
             })
        .def("get_pressures",
             [](const sph::ParticleSystem& ps) {
                 return py::array_t<float>(ps.size(), ps.get_pressures().data());
             });

    // Performance Stats
    py::class_<sph::SPHEngine::PerformanceStats>(m, "PerformanceStats")
        .def_readonly("total_time", &sph::SPHEngine::PerformanceStats::total_time)
        .def_readonly("neighbor_search_time", &sph::SPHEngine::PerformanceStats::neighbor_search_time)
        .def_readonly("density_computation_time", &sph::SPHEngine::PerformanceStats::density_computation_time)
        .def_readonly("force_computation_time", &sph::SPHEngine::PerformanceStats::force_computation_time)
        .def_readonly("integration_time", &sph::SPHEngine::PerformanceStats::integration_time)
        .def_readonly("max_neighbors", &sph::SPHEngine::PerformanceStats::max_neighbors)
        .def_readonly("total_neighbor_queries", &sph::SPHEngine::PerformanceStats::total_neighbor_queries);

    // Main SPH Engine
    py::class_<sph::SPHEngine>(m, "Simulator")
        .def(py::init<size_t>(), py::arg("max_particles") = 1000000)
        .def("initialize", &sph::SPHEngine::initialize, py::arg("params"))
        .def("initialize_dam_break", &sph::SPHEngine::initialize_dam_break)
        .def("initialize_fluid_drop", &sph::SPHEngine::initialize_fluid_drop)
        .def("initialize_granular_flow", &sph::SPHEngine::initialize_granular_flow)
        .def("add_particles", &sph::SPHEngine::add_particles)
        .def("clear_particles", &sph::SPHEngine::clear_particles)
        .def("step", &sph::SPHEngine::step, py::arg("dt") = 0.0f)
        .def("run_steps", &sph::SPHEngine::run_steps, py::arg("num_steps"), py::arg("adaptive_timestep") = true)
        .def("get_particles", &sph::SPHEngine::get_particles, py::return_value_policy::reference)
        .def("get_parameters", &sph::SPHEngine::get_parameters)
        .def("get_current_time", &sph::SPHEngine::get_current_time)
        .def("get_step_count", &sph::SPHEngine::get_step_count)
        .def("set_parameters", &sph::SPHEngine::set_parameters)
        .def("set_gravity", &sph::SPHEngine::set_gravity)
        .def("set_viscosity", &sph::SPHEngine::set_viscosity)
        .def("set_smoothing_length", &sph::SPHEngine::set_smoothing_length)
        .def("set_boundaries", &sph::SPHEngine::set_boundaries)
        .def("get_performance_stats", &sph::SPHEngine::get_performance_stats)
        .def("reset_performance_stats", &sph::SPHEngine::reset_performance_stats)
        .def("compute_conservation_errors",
             [](sph::SPHEngine& engine) {
                 float mass_error, energy_error;
                 engine.compute_conservation_errors(mass_error, energy_error);
                 return py::make_tuple(mass_error, energy_error);
             })
        .def("get_total_mass", &sph::SPHEngine::get_total_mass)
        .def("get_total_energy", &sph::SPHEngine::get_total_energy)
        .def("validate_simulation", &sph::SPHEngine::validate_simulation)
        .def("is_initialized", &sph::SPHEngine::is_initialized)

        // Numpy array accessors
        .def("get_positions",
             [](const sph::SPHEngine& engine) {
                 const auto& positions = engine.get_positions();
                 return py::array_t<float>({positions.size(), 3UL},
                     {3 * sizeof(float), sizeof(float)},
                     reinterpret_cast<const float*>(positions.data()));
             })
        .def("get_velocities",
             [](const sph::SPHEngine& engine) {
                 const auto& velocities = engine.get_velocities();
                 return py::array_t<float>({velocities.size(), 3UL},
                     {3 * sizeof(float), sizeof(float)},
                     reinterpret_cast<const float*>(velocities.data()));
             })
        .def("get_densities",
             [](const sph::SPHEngine& engine) {
                 const auto& densities = engine.get_densities();
                 return py::array_t<float>(densities.size(), densities.data());
             })
        .def("get_pressures",
             [](const sph::SPHEngine& engine) {
                 const auto& pressures = engine.get_pressures();
                 return py::array_t<float>(pressures.size(), pressures.data());
             });

    // Utility functions
    m.def("create_fluid_block",
          [](const py::array_t<float>& center, const py::array_t<float>& size, float spacing, float mass) {
              if (center.size() != 3 || size.size() != 3) {
                  throw std::invalid_argument("center and size must be 3D vectors");
              }
              glm::vec3 c(center.at(0), center.at(1), center.at(2));
              glm::vec3 s(size.at(0), size.at(1), size.at(2));
              return sph::create_fluid_block(c, s, spacing, mass);
          },
          py::arg("center"), py::arg("size"), py::arg("spacing"), py::arg("mass") = 1.0f);

    m.def("create_boundary_box",
          [](const py::array_t<float>& center, const py::array_t<float>& size, float spacing, float mass) {
              if (center.size() != 3 || size.size() != 3) {
                  throw std::invalid_argument("center and size must be 3D vectors");
              }
              glm::vec3 c(center.at(0), center.at(1), center.at(2));
              glm::vec3 s(size.at(0), size.at(1), size.at(2));
              return sph::create_boundary_box(c, s, spacing, mass);
          },
          py::arg("center"), py::arg("size"), py::arg("spacing"), py::arg("mass") = 1.0f);

    // Version info
    m.attr("__version__") = "1.0.0";
}

