# Breeze.jl

Documentation for Breeze.jl

## Overview

Breeze.jl is a Julia package for atmospheric thermodynamics calculations, particularly focused on Large Eddy Simulation (LES) applications.

## Features

- Accurate thermodynamic calculations
- Saturation vapor pressure computation
- Specific humidity calculations

## Installation

To install Breeze.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("Breeze")
```

## Quick Start

Here's a basic example of computing saturation specific humidity:

```julia
using Breeze

# Create thermodynamics instance
thermo = Breeze.AtmosphereThermodynamics()

# Calculate saturation specific humidity
T = 293.15  # Temperature in Kelvin (20°C)
ρ = 1.2     # Density in kg/m³
q★ = Breeze.saturation_specific_humidity(T, ρ, thermo, thermo.condensation)
``` 