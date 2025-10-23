#!/usr/bin/env python3
"""
SPH Particle Simulator - Python Visualization

This script demonstrates how to use the SPH simulator from Python
for scientific computing and visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import sph  # Our C++ module
import time

class SPHVisualizer:
    def __init__(self, simulator):
        self.simulator = simulator
        self.fig = plt.figure(figsize=(12, 8))

        # Create subplots
        self.ax1 = self.fig.add_subplot(221, projection='3d')  # 3D particle view
        self.ax2 = self.fig.add_subplot(222)  # Velocity histogram
        self.ax3 = self.fig.add_subplot(223)  # Density plot
        self.ax4 = self.fig.add_subplot(224)  # Energy conservation

        # Data storage for time series
        self.times = []
        self.energies = []
        self.mass_errors = []

        plt.tight_layout()

    def update_plot(self, frame):
        """Update all subplots for animation"""
        # Run simulation step
        self.simulator.step()

        # Get current data
        positions = np.array(self.simulator.get_positions())
        velocities = np.array(self.simulator.get_velocities())
        densities = np.array(self.simulator.get_densities())

        current_time = self.simulator.get_current_time()

        # Update time series data
        self.times.append(current_time)
        energy = self.simulator.get_total_energy()
        self.energies.append(energy)

        mass_error, _ = self.simulator.compute_conservation_errors()
        self.mass_errors.append(mass_error * 100.0)  # Convert to percentage

        # Clear all axes
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()

        # 3D Particle visualization
        if len(positions) > 0:
            speed = np.linalg.norm(velocities, axis=1)
            scatter = self.ax1.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                                     c=speed, cmap='plasma', s=2, alpha=0.8)
            self.ax1.set_xlabel('X')
            self.ax1.set_ylabel('Y')
            self.ax1.set_zlabel('Z')
            self.ax1.set_title(f'Particle Positions (t={current_time:.3f})')
            self.ax1.set_xlim([-1, 1])
            self.ax1.set_ylim([-1, 1])
            self.ax1.set_zlim([-1, 1])
            plt.colorbar(scatter, ax=self.ax1, shrink=0.8, label='Speed')

        # Velocity histogram
        if len(velocities) > 0:
            speed = np.linalg.norm(velocities, axis=1)
            self.ax2.hist(speed, bins=30, alpha=0.7, color='blue', edgecolor='black')
            self.ax2.set_xlabel('Speed')
            self.ax2.set_ylabel('Count')
            self.ax2.set_title('Velocity Distribution')
            self.ax2.grid(True, alpha=0.3)

        # Density plot (2D projection)
        if len(positions) > 0 and len(densities) > 0:
            self.ax3.scatter(positions[:, 0], positions[:, 1], c=densities,
                           cmap='viridis', s=1, alpha=0.6)
            self.ax3.set_xlabel('X')
            self.ax3.set_ylabel('Y')
            self.ax3.set_title('Density Distribution')
            self.ax3.set_xlim([-1, 1])
            self.ax3.set_ylim([-1, 1])
            plt.colorbar(self.ax3.collections[0], ax=self.ax3, shrink=0.8, label='Density')

        # Energy conservation plot
        if len(self.times) > 1:
            self.ax4.plot(self.times, self.energies, 'b-', label='Total Energy', linewidth=2)
            self.ax4.set_xlabel('Time')
            self.ax4.set_ylabel('Energy')
            self.ax4.set_title('Energy Conservation')
            self.ax4.grid(True, alpha=0.3)
            self.ax4.legend()

            # Add mass conservation as text
            if self.mass_errors:
                mass_err = self.mass_errors[-1]
                self.ax4.text(0.02, 0.98, f'Mass Error: {mass_err:.4f}%',
                            transform=self.ax4.transAxes, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        plt.tight_layout()
        return self.ax1, self.ax2, self.ax3, self.ax4

def run_dam_break_simulation(num_particles=5000, duration=2.0):
    """Run a dam break simulation with visualization"""
    print("Initializing SPH dam break simulation...")

    # Create simulator
    simulator = sph.Simulator(num_particles * 2)  # Extra space for boundaries

    # Configure parameters
    params = sph.SPHParameters()
    params.rest_density = 1000.0
    params.gas_constant = 2000.0
    params.viscosity = 0.001
    params.smoothing_length = 0.02
    params.particle_mass = 0.001
    params.timestep = 0.001
    params.gravity = -9.81

    simulator.initialize(params)
    simulator.initialize_dam_break()

    print(f"Initialized with {simulator.get_particles().size()} particles")

    # Create visualizer
    visualizer = SPHVisualizer(simulator)

    # Calculate number of frames
    dt = params.timestep
    num_frames = int(duration / dt)

    # Create animation
    ani = animation.FuncAnimation(visualizer.fig, visualizer.update_plot,
                                frames=num_frames, interval=50, blit=False)

    plt.show()

    # Print final statistics
    stats = simulator.get_performance_stats()
    print("
Simulation Statistics:")
    print(".3f")
    print(".1f")
    print(f"  Neighbor searches: {stats.total_neighbor_queries}")
    print(f"  Max neighbors per particle: {stats.max_neighbors}")

def run_performance_test(particle_counts=[1000, 5000, 10000, 25000]):
    """Run performance tests for different particle counts"""
    print("Running SPH performance tests...")

    results = []

    for num_particles in particle_counts:
        print(f"\nTesting with {num_particles} particles...")

        # Create simulator
        simulator = sph.Simulator(num_particles * 2)

        # Configure for performance
        params = sph.SPHParameters()
        params.smoothing_length = 0.02
        params.particle_mass = 0.001

        simulator.initialize(params)
        simulator.initialize_dam_break()

        # Warm up
        for _ in range(10):
            simulator.step()

        # Time measurement
        start_time = time.time()
        num_steps = 100

        for _ in range(num_steps):
            simulator.step()

        end_time = time.time()
        total_time = end_time - start_time

        fps = num_steps / total_time
        ms_per_frame = (total_time / num_steps) * 1000.0

        results.append({
            'particles': num_particles,
            'fps': fps,
            'ms_per_frame': ms_per_frame
        })

        print(".2f")
        print(".2f")

    # Plot results
    particles = [r['particles'] for r in results]
    fps_values = [r['fps'] for r in results]

    plt.figure(figsize=(10, 6))
    plt.loglog(particles, fps_values, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('Number of Particles')
    plt.ylabel('FPS')
    plt.title('SPH Performance Scaling')
    plt.grid(True, alpha=0.3)

    # Add trend line
    log_particles = np.log(particles)
    log_fps = np.log(fps_values)
    coeffs = np.polyfit(log_particles, log_fps, 1)
    trend_particles = np.logspace(np.log10(min(particles)), np.log10(max(particles)), 100)
    trend_fps = np.exp(coeffs[1]) * np.power(trend_particles, coeffs[0])
    plt.loglog(trend_particles, trend_fps, 'r--', alpha=0.7, label=f'Slope: {coeffs[0]:.2f}')

    plt.legend()
    plt.tight_layout()
    plt.show()

def create_particle_data_export():
    """Demonstrate data export for analysis"""
    print("Creating particle data export example...")

    simulator = sph.Simulator(10000)
    params = sph.SPHParameters()
    simulator.initialize(params)
    simulator.initialize_dam_break()

    # Run a few steps
    for _ in range(50):
        simulator.step()

    # Export data
    positions = np.array(simulator.get_positions())
    velocities = np.array(simulator.get_velocities())
    densities = np.array(simulator.get_densities())
    pressures = np.array(simulator.get_pressures())

    print("Exported data shapes:")
    print(f"  Positions: {positions.shape}")
    print(f"  Velocities: {velocities.shape}")
    print(f"  Densities: {densities.shape}")
    print(f"  Pressures: {pressures.shape}")

    # Save to numpy format
    np.savez('sph_data_export.npz',
             positions=positions,
             velocities=velocities,
             densities=densities,
             pressures=pressures,
             time=simulator.get_current_time())

    print("Data saved to 'sph_data_export.npz'")

    # Example analysis
    speed = np.linalg.norm(velocities, axis=1)
    print("
Analysis:")
    print(".3f")
    print(".1f")
    print(".1f")
    print(".1f")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SPH Particle Simulator - Python Interface')
    parser.add_argument('--mode', choices=['visualize', 'performance', 'export'],
                       default='visualize', help='Operation mode')
    parser.add_argument('--particles', type=int, default=5000,
                       help='Number of particles for simulation')
    parser.add_argument('--duration', type=float, default=2.0,
                       help='Simulation duration in seconds')

    args = parser.parse_args()

    if args.mode == 'visualize':
        run_dam_break_simulation(args.particles, args.duration)
    elif args.mode == 'performance':
        run_performance_test()
    elif args.mode == 'export':
        create_particle_data_export()

