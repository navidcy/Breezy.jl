module AquaSkyLES

export MoistAirBuoyancy

using Oceananigans

import Oceananigans.BuoyancyFormulations: AbstractBuoyancyFormulation,
                                          buoyancy_perturbationᶜᶜᶜ,
                                          required_tracers


struct MoistAirBuoyancy{FT} <: AbstractBuoyancyFormulation{Nothing}
    expansion_coefficient :: FT
    reference_potential_temperature :: FT
    gas_constant_ratio :: FT
end

function MoistAirBuoyancy(FT=Oceananigans.defaults.FloatType;
                          expansion_coefficient = 3.27e-2,
                          reference_potential_temperature = 0,
                          gas_constant_ratio = 1.61)

    return MoistAirBuoyancy{FT}(expansion_coefficient,
                                reference_potential_temperature,
                                gas_constant_ratio)
end

required_tracers(::MoistAirBuoyancy) = (:θ, :q)

@inline function buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, mb::MoistAirBuoyancy, tracers)
    β = mb.expansion_coefficient
    θ₀ = mb.reference_potential_temperature
    ϵᵥ = mb.gas_constant_ratio
    δ = ϵᵥ - 1
    θ = @inbounds tracers.θ[i, j, k]
    q = @inbounds tracers.q[i, j, k]
    θᵥ = θ * (1 + δ * q)
    return β * (θᵥ - θ₀)
end

end # module AquaSkyLES
