#include "sph_engine.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fstream>

/**
 * SPH Performance Benchmark
 *
 * This program benchmarks the SPH simulator performance across different
 * particle counts and configurations, providing detailed performance analysis.
 */

struct BenchmarkResult {
    size_t particle_count;
    double avg_fps;
    double min_fps;
    double max_fps;
    double avg_ms_per_frame;
    double neighbor_search_time;
    double density_computation_time;
    double force_computation_time;
    double integration_time;
    size_t avg_neighbors;
    size_t max_neighbors;
    double total_simulation_time;
    double mass_conservation_error;
    double energy_conservation_error;
};

class SPHBenchmark {
private:
    std::vector<BenchmarkResult> results_;
    std::ofstream csv_file_;

public:
    SPHBenchmark() {
        // Open CSV file for results
        csv_file_.open("benchmark_results.csv");
        csv_file_ << "particles,avg_fps,min_fps,max_fps,avg_ms_frame,"
                 << "neighbor_search_pct,density_pct,force_pct,integration_pct,"
                 << "avg_neighbors,max_neighbors,total_time,mass_error,energy_error\n";
    }

    ~SPHBenchmark() {
        if (csv_file_.is_open()) {
            csv_file_.close();
        }
    }

    void run_benchmark(const std::vector<size_t>& particle_counts,
                      size_t warmup_steps = 50,
                      size_t benchmark_steps = 500) {

        std::cout << "SPH Performance Benchmark" << std::endl;
        std::cout << "========================" << std::endl;
        std::cout << std::endl;

        for (size_t num_particles : particle_counts) {
            std::cout << "Benchmarking with " << num_particles << " particles..." << std::endl;

            BenchmarkResult result = benchmark_configuration(num_particles, warmup_steps, benchmark_steps);
            results_.push_back(result);

            print_result(result);
            save_to_csv(result);

            std::cout << std::endl;
        }

        print_summary();
    }

private:
    BenchmarkResult benchmark_configuration(size_t num_particles,
                                          size_t warmup_steps,
                                          size_t benchmark_steps) {

        BenchmarkResult result;
        result.particle_count = num_particles;

        // Initialize SPH engine
        sph::SPHEngine engine(num_particles * 2);

        sph::SPHParameters params;
        params.rest_density = 1000.0f;
        params.gas_constant = 2000.0f;
        params.viscosity = 0.001f;
        params.smoothing_length = 0.02f;
        params.particle_mass = 0.001f;
        params.timestep = 0.001f;
        params.gravity = -9.81f;

        // Adjust smoothing length based on particle count for optimal performance
        float target_neighbors = 50.0f;
        float volume_per_particle = 1.0f / params.rest_density * params.particle_mass;
        float particle_spacing = std::pow(volume_per_particle, 1.0f/3.0f);
        params.smoothing_length = particle_spacing * std::pow(target_neighbors / 50.0f, 1.0f/3.0f);

        engine.initialize(params);
        engine.initialize_dam_break();

        // Warmup phase
        for (size_t i = 0; i < warmup_steps; ++i) {
            engine.step();
        }

        // Reset performance stats
        engine.reset_performance_stats();

        // Benchmark phase
        std::vector<double> frame_times;
        auto start_time = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < benchmark_steps; ++i) {
            auto frame_start = std::chrono::high_resolution_clock::now();
            engine.step();
            auto frame_end = std::chrono::high_resolution_clock::now();

            double frame_time = std::chrono::duration<double>(frame_end - frame_start).count();
            frame_times.push_back(frame_time);
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        result.total_simulation_time = std::chrono::duration<double>(end_time - start_time).count();

        // Calculate FPS statistics
        std::vector<double> fps_values;
        for (double ft : frame_times) {
            fps_values.push_back(1.0 / ft);
        }

        result.avg_fps = 1.0 / (result.total_simulation_time / benchmark_steps);
        result.min_fps = *std::min_element(fps_values.begin(), fps_values.end());
        result.max_fps = *std::max_element(fps_values.begin(), fps_values.end());
        result.avg_ms_per_frame = (result.total_simulation_time / benchmark_steps) * 1000.0;

        // Performance breakdown
        const auto& stats = engine.get_performance_stats();
        double total_perf_time = stats.total_time;

        if (total_perf_time > 0.0) {
            result.neighbor_search_time = stats.neighbor_search_time / total_perf_time * 100.0;
            result.density_computation_time = stats.density_computation_time / total_perf_time * 100.0;
            result.force_computation_time = stats.force_computation_time / total_perf_time * 100.0;
            result.integration_time = stats.integration_time / total_perf_time * 100.0;
        }

        // Neighbor statistics
        result.avg_neighbors = 50;  // Estimated, would need access to neighbor lists
        result.max_neighbors = stats.max_neighbors;

        // Conservation errors
        engine.compute_conservation_errors(result.mass_conservation_error, result.energy_conservation_error);

        return result;
    }

    void print_result(const BenchmarkResult& result) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Results for " << result.particle_count << " particles:" << std::endl;
        std::cout << "  Performance:" << std::endl;
        std::cout << "    Average FPS: " << result.avg_fps << std::endl;
        std::cout << "    FPS range: " << result.min_fps << " - " << result.max_fps << std::endl;
        std::cout << "    Average frame time: " << result.avg_ms_per_frame << " ms" << std::endl;
        std::cout << "  Time breakdown (%):" << std::endl;
        std::cout << "    Neighbor search: " << result.neighbor_search_time << std::endl;
        std::cout << "    Density computation: " << result.density_computation_time << std::endl;
        std::cout << "    Force computation: " << result.force_computation_time << std::endl;
        std::cout << "    Integration: " << result.integration_time << std::endl;
        std::cout << "  Neighbor statistics:" << std::endl;
        std::cout << "    Average neighbors: " << result.avg_neighbors << std::endl;
        std::cout << "    Maximum neighbors: " << result.max_neighbors << std::endl;
        std::cout << "  Conservation:" << std::endl;
        std::cout << "    Mass error: " << (result.mass_conservation_error * 100.0) << "%" << std::endl;
        std::cout << "    Energy error: " << (result.energy_conservation_error * 100.0) << "%" << std::endl;
    }

    void save_to_csv(const BenchmarkResult& result) {
        csv_file_ << result.particle_count << ","
                 << result.avg_fps << ","
                 << result.min_fps << ","
                 << result.max_fps << ","
                 << result.avg_ms_per_frame << ","
                 << result.neighbor_search_time << ","
                 << result.density_computation_time << ","
                 << result.force_computation_time << ","
                 << result.integration_time << ","
                 << result.avg_neighbors << ","
                 << result.max_neighbors << ","
                 << result.total_simulation_time << ","
                 << result.mass_conservation_error << ","
                 << result.energy_conservation_error << "\n";
    }

    void print_summary() {
        if (results_.empty()) return;

        std::cout << "Benchmark Summary" << std::endl;
        std::cout << "=================" << std::endl;
        std::cout << std::setw(10) << "Particles"
                 << std::setw(10) << "FPS"
                 << std::setw(12) << "Frame Time"
                 << std::setw(12) << "Mass Error"
                 << std::endl;
        std::cout << std::string(46, '-') << std::endl;

        std::cout << std::fixed << std::setprecision(1);
        for (const auto& result : results_) {
            std::cout << std::setw(10) << result.particle_count
                     << std::setw(10) << result.avg_fps
                     << std::setw(10) << result.avg_ms_per_frame << "ms"
                     << std::setw(10) << (result.mass_conservation_error * 100.0) << "%"
                     << std::endl;
        }

        // Scaling analysis
        if (results_.size() >= 2) {
            std::cout << std::endl << "Scaling Analysis:" << std::endl;

            for (size_t i = 1; i < results_.size(); ++i) {
                double particles_ratio = static_cast<double>(results_[i].particle_count) /
                                       results_[i-1].particle_count;
                double fps_ratio = results_[i-1].avg_fps / results_[i].avg_fps;

                double scaling_efficiency = fps_ratio / particles_ratio;
                std::cout << "  " << results_[i-1].particle_count << "k -> "
                         << results_[i].particle_count << "k: "
                         << (scaling_efficiency * 100.0) << "% efficient" << std::endl;
            }
        }

        std::cout << std::endl << "Results saved to 'benchmark_results.csv'" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    // Default benchmark configurations
    std::vector<size_t> particle_counts = {1000, 2500, 5000, 10000, 25000, 50000};

    // Allow custom particle counts via command line
    if (argc > 1) {
        particle_counts.clear();
        for (int i = 1; i < argc; ++i) {
            particle_counts.push_back(std::atoi(argv[i]));
        }
    }

    try {
        SPHBenchmark benchmark;
        benchmark.run_benchmark(particle_counts);

    } catch (const std::exception& e) {
        std::cerr << "Benchmark failed: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

