module AquaSkyLES

export UnsaturatedMoistAirBuoyancy

using Oceananigans
using Oceananigans: AbstractModel
using Oceananigans.Grids: AbstractGrid

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
    reference_state :: ReferenceState{FT}
end

function MoistAirBuoyancy(FT=Oceananigans.defaults.FloatType;
                           thermodynamics = AtmosphereThermodynamics(FT),
                           reference_state = ReferenceState{FT}(101325, 290))

    return MoistAirBuoyancy{FT}(thermodynamics, reference_state)
end

required_tracers(::MoistAirBuoyancy) = (:θ, :q)

reference_density(z, mb::MoistAirBuoyancy) = reference_density(z, mb.reference_state, mb.thermodynamics)
base_density(mb::MoistAirBuoyancy) = base_density(mb.reference_state, mb.thermodynamics)

#####
##### 
#####

const c = Center()

@inline function buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, mb::MoistAirBuoyancy, tracers)
    z = Oceananigans.Grids.znode(i, j, k, grid, c, c, c)
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

const c = Center()

#####
##### Temperature
#####

function temperature(i, j, k, grid::AbstractGrid, mb::MoistAirBuoyancy, θ, q)
    z = Oceananigans.Grids.znode(i, j, k, grid, c, c, c)
    θi = @inbounds θ[i, j, k]
    qi = @inbounds q[i, j, k]
    𝒰 = ThermodynamicState(θi, qi, z)
    return temperature(𝒰, mb.reference_state, mb.thermodynamics)
end

struct TemperatureKernelFunction end

@inline (::TemperatureKernelFunction)(i, j, k, grid, buoyancy, θ, q) =
    temperature(i, j, k, grid, buoyancy, θ, q)

function TemperatureField(model)
    func = TemperatureKernelFunction()
    grid = model.grid
    buoyancy = model.buoyancy.formulation
    θ = model.tracers.θ
    q = model.tracers.q
    op = KernelFunctionOperation{Center, Center, Center}(func, grid, buoyancy, θ, q)
    return Field(op)
end

#####
##### Saturation specific humidity
#####

@inline function saturation_specific_humidity(i, j, k, grid, mb::MoistAirBuoyancy, T)
    z = Oceananigans.Grids.znode(i, j, k, grid, c, c, c)
    Ti = @inbounds T[i, j, k]
    return saturation_specific_humidity(Ti, z, mb.reference_state, mb.thermodynamics)
end

struct SaturationKernelFunction{T}
    temperature :: T
end

@inline function (kernel::SaturationKernelFunction)(i, j, k, grid, buoyancy)
    T = kernel.temperature
    return saturation_specific_humidity(i, j, k, grid, buoyancy, T)
end

SaturationKernelFunction(model::AbstractModel, T=TemperatureField(model)) =
    SaturationKernelFunction(T)

function SaturationField(model, T=TemperatureField(model))
    func = SaturationKernelFunction(model, T)
    grid = model.grid
    buoyancy = model.buoyancy.formulation
    op = KernelFunctionOperation{Center, Center, Center}(func, grid, buoyancy)
    return Field(op)
end

#####
##### Condensate
#####

struct CondensateKernelFunction{T}
    temperature :: T
end

@inline function condensate_specific_humidity(i, j, k, grid, mb::MoistAirBuoyancy, T, q)
    z = Oceananigans.Grids.znode(i, j, k, grid, c, c, c)
    Ti = @inbounds T[i, j, k]
    qi = @inbounds q[i, j, k]
    qˡ = condensate_specific_humidity(Ti, qi, z, mb.reference_state, mb.thermodynamics)
    return qˡ
end

@inline function (kernel::CondensateKernelFunction)(i, j, k, grid, buoyancy, q)
    T = kernel.temperature
    return condensate_specific_humidity(i, j, k, grid, buoyancy, T, q)
end

CondensateKernelFunction(model::AbstractModel, T=TemperatureField(model)) =
    CondensateKernelFunction(T)

function CondensateField(model, T=TemperatureField(model))
    func = CondensateKernelFunction(model, T)
    grid = model.grid
    buoyancy = model.buoyancy.formulation
    q = model.tracers.q
    op = KernelFunctionOperation{Center, Center, Center}(func, grid, buoyancy, q)
    return Field(op)
end

end # module AquaSkyLES
