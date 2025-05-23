using Oceananigans.Utils: prettysummary

using ..Thermodynamics:
    AtmosphereThermodynamics,
    ReferenceConstants,
    reference_pressure,
    reference_density,
    mixture_gas_constant,
    mixture_heat_capacity,
    dry_air_gas_constant

struct AnelasticFormulation{FT, F}
    constants :: ReferenceConstants{FT}
    reference_pressure :: F
    reference_density :: F
end

function Base.summary(formulation::AnelasticFormulation)
    p₀ = formulation.constants.base_pressure
    θᵣ = formulation.constants.reference_potential_temperature
    return string("AnelasticFormulation(p₀=", prettysummary(p₀), 
                  ", θᵣ=", prettysummary(θᵣ), ")")
end

Base.show(io::IO, formulation::AnelasticFormulation) = print(io, "AnelasticFormulation")

field_names(::AnelasticFormulation, tracer_names) = (:ρu, :ρv, :ρw, :e, :ρq, tracer_names...)

struct AnelasticThermodynamicState{FT}
    potential_temperature :: FT
    specific_humidity :: FT
    reference_density :: FT
    reference_pressure :: FT
    exner_function :: FT
end

function AnelasticFormulation(grid, state_constants, thermo)
    pᵣ = Field{Nothing, Nothing, Center}(grid)
    ρᵣ = Field{Nothing, Nothing, Center}(grid)
    set!(pᵣ, z -> reference_pressure(z, state_constants, thermo))
    set!(ρᵣ, z -> reference_density(z, state_constants, thermo))
    fill_halo_regions!(pᵣ)
    fill_halo_regions!(ρᵣ)
    return AnelasticFormulation(state_constants, pᵣ, ρᵣ)
end

function thermodynamic_state(i, j, k, grid, formulation::AnelasticFormulation, thermo, energy, absolute_humidity)
    @inbounds begin
        e = energy[i, j, k]
        pᵣ = formulation.reference_pressure[i, j, k]
        ρᵣ = formulation.reference_density[i, j, k]
        ρq = absolute_humidity[i, j, k]
    end

    cᵖᵈ = thermo.dry_air.heat_capacity
    θ = e / (cᵖᵈ * ρᵣ)

    q = ρq / ρᵣ
    Rᵐ = mixture_gas_constant(q, thermo)
    cᵖᵐ = mixture_heat_capacity(q, thermo)

    p₀ = formulation.constants.base_pressure
    Π = (pᵣ / p₀)^(Rᵐ / cᵖᵐ)

    return AnelasticThermodynamicState(θ, q, ρᵣ, pᵣ, Π)
end

@inline function specific_volume(i, j, k, grid, formulation, temperature, specific_humidity, thermo)
    @inbounds begin
        q  = specific_humidity[i, j, k]
        pᵣ = formulation.reference_pressure[i, j, k]
        T = temperature[i, j, k]
    end

    Rᵐ = mixture_gas_constant(q, thermo)

    return Rᵐ * T / pᵣ
end

@inline function reference_specific_volume(i, j, k, grid, formulation, thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    pᵣ = @inbounds formulation.reference_pressure[i, j, k]
    θᵣ = formulation.constants.reference_potential_temperature
    return Rᵈ * θᵣ / pᵣ
end

function collect_prognostic_fields(::AnelasticFormulation,
                                   density,
                                   momentum,
                                   energy,
                                   absolute_humidity,
                                   condensates,
                                   tracers)

    thermodynamic_variables = (e=energy, ρq=absolute_humidity)

    return merge(momentum, thermodynamic_variables, condensates, tracers)
end

function materialize_momentum_and_velocities(formulation::AnelasticFormulation, grid, boundary_conditions)
    ρu = XFaceField(grid, boundary_conditions=boundary_conditions.ρu)
    ρv = YFaceField(grid, boundary_conditions=boundary_conditions.ρv)
    ρw = ZFaceField(grid, boundary_conditions=boundary_conditions.ρw)
    momentum = (; ρu, ρv, ρw)

    velocity_bcs = NamedTuple(name => FieldBoundaryConditions() for name in (:u, :v, :w))
    velocity_bcs = regularize_field_boundary_conditions(velocity_bcs, grid, (:u, :v, :w))
    u = XFaceField(grid, boundary_conditions=velocity_bcs.u)
    v = YFaceField(grid, boundary_conditions=velocity_bcs.v)
    w = ZFaceField(grid, boundary_conditions=velocity_bcs.w)
    velocities = (; u, v, w)

    return velocities, momentum
end
