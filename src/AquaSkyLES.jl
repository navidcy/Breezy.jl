module AquaSkyLES

export UnsaturatedMoistAirBuoyancy

using Oceananigans

import Oceananigans.BuoyancyFormulations: AbstractBuoyancyFormulation,
                                          buoyancy_perturbationᶜᶜᶜ,
                                          required_tracers

include("atmospheric_thermodynamics.jl")

struct UnsaturatedMoistAirBuoyancy{FT} <: AbstractBuoyancyFormulation{Nothing}
    expansion_coefficient :: FT
    reference_potential_temperature :: FT
    gas_constant_ratio :: FT
end

function UnsaturatedMoistAirBuoyancy(FT=Oceananigans.defaults.FloatType;
                                     expansion_coefficient = 3.27e-2,
                                     reference_potential_temperature = 0,
                                     gas_constant_ratio = 1.61)

    return UnsaturatedMoistAirBuoyancy{FT}(expansion_coefficient,
                                           reference_potential_temperature,
                                           gas_constant_ratio)
end

required_tracers(::UnsaturatedMoistAirBuoyancy) = (:θ, :q)

@inline function buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, mb::UnsaturatedMoistAirBuoyancy, tracers)
    β = mb.expansion_coefficient
    θ₀ = mb.reference_potential_temperature
    ϵᵥ = mb.gas_constant_ratio
    δ = ϵᵥ - 1
    θ = @inbounds tracers.θ[i, j, k]
    q = @inbounds tracers.q[i, j, k]
    θᵥ = θ * (1 + δ * q)
    return β * (θᵥ - θ₀)
end


struct MoistAirBuoyancy{FT} <: AbstractBuoyancyFormulation{Nothing}
    thermodynamics :: AtmosphereThermodynamics{FT}
    reference_state :: FT
end

function MoistAirBuoyancy(FT=Oceananigans.defaults.FloatType;
                           thermodynamics = AtmosphereThermodynamics(FT),
                           reference_state = ReferenceState{FT}(101325, 20))

    return MoistAirBuoyancy{FT}(thermodynamics, reference_state)
end

required_tracers(::MoistAirBuoyancy) = (:θ, :q)

const c = Center()

@inline function buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, mb::MoistAirBuoyancy, tracers)
    z = znode(i, j, k, grid, c, c, c)
    θ = @inbounds tracers.θ[i, j, k]
    q = @inbounds tracers.q[i, j, k]
    𝒰 = ThermodynamicState(θ, q, z)

    ρ₀ = base_density(mb.reference_state, mb.thermodynamics)
    αʳ = reference_specific_volume(z, mb.reference_state, mb.thermodynamics)
    g = mb.thermodynamics.gravitational_acceleration

    # Perform saturation adjustment
    α = specific_volume(𝒰, mb.reference_state, mb.thermodynamics)

    return ρ₀ * g * (α - αʳ)
end

end # module AquaSkyLES
