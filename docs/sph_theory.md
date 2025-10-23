# SPH Theory and Implementation

## Overview

Smoothed Particle Hydrodynamics (SPH) is a computational method used for simulating fluid flows and continuum mechanics problems. Unlike traditional grid-based methods (FEM, FVM), SPH uses a Lagrangian approach where the computational domain is discretized into particles that carry physical properties and move with the material.

## Fundamental Principles

### Kernel Approximation

SPH approximates any field quantity A(r) using:

```
A(r) = ∫ A(r') W(r - r', h) dr'
```

Where:
- `W(r, h)` is the smoothing kernel
- `h` is the smoothing length
- The integral is approximated by summation over neighboring particles

### Particle Approximation

For a system of N particles, the integral becomes:

```
Aᵢ = ∑ⱼ mⱼ Aⱼ / ρⱼ W(rᵢⱼ, h)
```

Where:
- `mⱼ` is particle mass
- `ρⱼ` is particle density
- `rᵢⱼ = rᵢ - rⱼ` is the relative position vector

## SPH Equations

### Density Computation

The density at particle i is:

```
ρᵢ = ∑ⱼ mⱼ W(rᵢⱼ, h)
```

### Pressure Computation

Using the ideal gas equation of state:

```
Pᵢ = k(ρᵢ - ρ₀)
```

Where:
- `k` is the gas constant
- `ρ₀` is the rest density

### Momentum Equation (Navier-Stokes)

The acceleration of particle i includes pressure, viscosity, and external forces:

```
d²rᵢ/dt² = -∑ⱼ mⱼ (Pᵢ + Pⱼ)/(2ρⱼ) ∇W(rᵢⱼ, h) + μ ∑ⱼ mⱼ (vⱼ - vᵢ)/ρⱼ ∇²W(rᵢⱼ, h) + g
```

#### Pressure Force

```
Fᵢᵖʳᵉˢˢᵘʳᵉ = -∑ⱼ mⱼ (Pᵢ + Pⱼ)/(2ρⱼ) ∇W(rᵢⱼ, h)
```

#### Viscosity Force

```
Fᵢᵛⁱˢᶜᵒˢⁱᵗʸ = μ ∑ⱼ mⱼ (vⱼ - vᵢ)/ρⱼ ∇²W(rᵢⱼ, h)
```

## Smoothing Kernels

### Cubic Spline Kernel (Most Common)

The cubic spline kernel is widely used due to its smoothness and compact support:

```
W(q) = σ * (2/3 - q² + 1/2 q³)    for 0 ≤ q < 1
W(q) = σ * (1/6)(2 - q)³         for 1 ≤ q < 2
W(q) = 0                         for q ≥ 2
```

Where:
- `q = |r| / h`
- `σ` is the normalization constant (1/(πh³) in 3D)

### Kernel Properties

A good smoothing kernel should satisfy:
1. **Normalization**: ∫ W(r,h) dr = 1
2. **Compact support**: W(r,h) = 0 for |r| ≥ 2h
3. **Monotonicity**: dW/dq < 0 for q > 0
4. **Continuity**: Continuous first and second derivatives

## Time Integration

### Verlet Integration (Recommended)

```
r(t+Δt) = 2r(t) - r(t-Δt) + a(t) Δt²
v(t+Δt) = [r(t+Δt) - r(t-Δt)] / (2Δt)
```

### Leapfrog Integration

```
v(t+½Δt) = v(t-½Δt) + a(t) Δt
r(t+Δt) = r(t) + v(t+½Δt) Δt
v(t+Δt) = v(t+½Δt) + a(t+Δt) Δt
```

## Boundary Conditions

### Repulsive Boundaries

Particles near boundaries experience repulsive forces:

```
Fᵇᵒᵘⁿᵈ = -k (1/r - 1/r₀) r̂    for r < r₀
```

### Periodic Boundaries

Particles crossing boundaries wrap around:

```
x = x - L * floor(x/L)
```

## Stability and Accuracy

### CFL Condition

The timestep must satisfy:

```
Δt ≤ 0.4 * h / max|v|
```

### Artificial Viscosity

Added to prevent particle penetration:

```
Πᵢⱼ = -α cᵢⱼ μᵢⱼ / ρ̄ᵢⱼ
```

Where:
- `μᵢⱼ = h (vᵢ - vⱼ) · (rᵢ - rⱼ) / (|rᵢ - rⱼ|² + η²)`
- `cᵢⱼ = (cᵢ + cⱼ)/2`
- `ρ̄ᵢⱼ = (ρᵢ + ρⱼ)/2`

## SPH Variants

### Weakly Compressible SPH (WCSPH)

Uses stiff equation of state for incompressible flows:

```
P = (ρ/ρ₀)^γ - 1
```

### Incompressible SPH (ISPH)

Solves pressure Poisson equation:

```
∇²P = ρ / Δt² ∑ⱼ ∇Wᵢⱼ · (vⱼ - vᵢ)
```

### Moving Least Squares SPH (MLSPH)

Uses higher-order interpolation for better accuracy.

## Applications

- Free surface flows
- Multi-phase flows
- Fluid-structure interaction
- Granular flows
- Astrophysics (star formation, galaxy dynamics)
- Environmental engineering
- Biomechanics

## References

1. Lucy, L. B. (1977). A numerical approach to the testing of the fission hypothesis.
2. Gingold, R. A., & Monaghan, J. J. (1977). Smoothed particle hydrodynamics.
3. Monaghan, J. J. (1992). Smoothed particle hydrodynamics.
4. Liu, G. R., & Liu, M. B. (2003). Smoothed Particle Hydrodynamics: A Meshfree Method.

