using ..MoistThermodynamics:
    AtmosphereThermodynamics,
    ReferenceConstants,
    reference_pressure,
    reference_density,
    mixture_gas_constant,
    mixture_heat_capacity,
    dry_air_gas_constant

using Oceananigans
using Oceananigans.Advection: Centered, adapt_advection_order
using Oceananigans.Grids: ZDirection
using Oceananigans.Models: AbstractModel
using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.TimeSteppers: TimeStepper
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, regularize_field_boundary_conditions
using Oceananigans.Solvers: FourierTridiagonalPoissonSolver
using Oceananigans.Utils: launch!

using KernelAbstractions: @kernel, @index

struct AnelasticReferenceState{FT, F}
    constants :: ReferenceConstants{FT}
    pressure :: F
    density :: F
end

struct AnelasticParcelState{FT}
    potential_temperature :: FT
    specific_humidity :: FT
    reference_pressure :: FT
    reference_density :: FT
    exner_function :: FT
end

@inline function AnelasticParcelState(θ, q, pʳ, ρʳ, thermo::AtmosphereThermodynamics)
    Rᵐ = mixture_gas_constant(q, thermo)
    cᵖᵐ = mixture_heat_capacity(q, thermo)
    inv_ϰᵐ = Rᵐ / cᵖᵐ
    Π = @inbounds (pʳ / pʳ[1])^inv_ϰᵐ
    return AnelasticParcelState(θ, q, pʳ, ρʳ, Π)
end

function AnelasticReferenceState(grid, ref_constants, thermo)
    pʳ = Field{Center, Nothing, Nothing}(grid)
    ρʳ = Field{Center, Nothing, Nothing}(grid)
    set!(pʳ, z -> reference_pressure(z, ref_constants, thermo))
    set!(ρʳ, z -> reference_density(z, ref_constants, thermo))
    return AnelasticReferenceState(ref_constants, pʳ, ρʳ)
end

@inline function anelastic_specific_volume(i, j, k, grid, reference_state, temperature, specific_humidity, thermo)
    @inbounds begin
        q  = specific_humidity[i, j, k]
        pʳ = reference_state.pressure[i, j, k]
        T = temperature[i, j, k]
    end

    Rᵐ = mixture_gas_constant(q, thermo)

    return Rᵐ * T / pʳ
end

@inline function reference_specific_volume(i, j, k, grid, reference_state, thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    pʳ = @inbounds reference_state.pressure[i, j, k]
    θʳ = reference_state.constants.reference_potential_temperature
    return Rᵈ * θʳ / pʳ
end

field_names(::AnelasticReferenceState, tracer_names) = (:ρu, :ρv, :ρw, :ρθ, :ρq, tracer_names...)

function collect_prognostic_fields(::AnelasticReferenceState,
                                   density,
                                   momentum,
                                   potential_temperature_density,
                                   absolute_humidity,
                                   condensates,
                                   tracers)

    thermodynamic_variables = (ρθ=potential_temperature_density, ρq=absolute_humidity)

    return merge(momentum, thermodynamic_variables, condensates, tracers)
end

function default_reference_state(grid, thermo)
    FT = eltype(grid)
    base_pressure = convert(FT, 101325)
    potential_temperature = convert(FT, 288)
    constants = ReferenceConstants(base_pressure, potential_temperature)
    return AnelasticReferenceState(grid, constants, thermo)
end

function materialize_momentum_and_velocities(reference_state::AnelasticReferenceState, grid, boundary_conditions)
    ρu = XFaceField(grid, boundary_conditions=boundary_conditions.ρu)
    ρv = YFaceField(grid, boundary_conditions=boundary_conditions.ρv)
    ρw = ZFaceField(grid, boundary_conditions=boundary_conditions.ρw)
    momentum = (; ρu, ρv, ρw)

    u = XFaceField(grid)
    v = YFaceField(grid)
    w = ZFaceField(grid)
    velocities = (; u, v, w)

    return velocities, momentum
end

materialize_condenstates(microphysics, grid) = NamedTuple() #(; qˡ=CenterField(grid), qᵛ=CenterField(grid))
materialize_density(reference_state, grid) = CenterField(grid)

struct WarmPhaseSaturationAdjustment end

tupleit(t::Tuple) = t
tupleit(t) = tuple(t)

mutable struct AtmosphereModel{Re, Ar, Ts, Gr, Cl, Th, De, Mo, Pd, Wa, En, Hu,
                               Te, Pr, Pp, So, Ve, Tr, Ad, Co, Mi, Cn} <: AbstractModel{Ts, Ar}
    architecture :: Ar
    grid :: Gr
    clock :: Cl
    reference_state :: Re
    thermodynamics :: Th
    density :: De
    momentum :: Mo
    potential_temperature_density :: Pd
    absolute_humidity :: Wa
    potential_temperature :: En # TODO: generalize energy formulation
    specific_humidity :: Hu
    temperature :: Te
    pressure :: Pr
    hydrostatic_pressure_anomaly :: Pp
    pressure_solver :: So
    velocities :: Ve
    tracers :: Tr
    advection :: Ad
    coriolis :: Co
    microphysics :: Mi
    condensates :: Cn
    timestepper :: Ts
end

function AtmosphereModel(grid;
                         clock = Clock(grid),
                         thermodynamics = AtmosphereThermodynamics(eltype(grid)),
                         reference_state = default_reference_state(grid, thermodynamics),
                         specific_humidity = CenterField(grid),
                         tracers = tuple(),
                         coriolis = nothing,
                         boundary_conditions = NamedTuple(),
                         advection = WENO(order=5),
                         microphysics = WarmPhaseSaturationAdjustment(),
                         timestepper = :RungeKutta3)

    arch = grid.architecture
    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)

    hydrostatic_pressure_anomaly = CenterField(grid)
    pressure = CenterField(grid)

    # Next, we form a list of default boundary conditions:
    names = field_names(reference_state, tracers)
    default_boundary_conditions = NamedTuple{names}(FieldBoundaryConditions() for _ in names)
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, names)

    density = materialize_density(reference_state, grid)
    velocities, momentum = materialize_momentum_and_velocities(reference_state, grid, boundary_conditions)
    tracers = NamedTuple(n => CenterField(grid, boundary_conditions=boundary_conditions[n]) for name in tracers)
    condensates = materialize_condenstates(microphysics, grid)
    advection = adapt_advection_order(advection, grid)

    potential_temperature_density = CenterField(grid, boundary_conditions=boundary_conditions.ρθ)
    absolute_humidity = CenterField(grid, boundary_conditions=boundary_conditions.ρq)
    potential_temperature = CenterField(grid)
    specific_humidity = CenterField(grid)
    temperature = CenterField(grid)

    prognostic_fields = collect_prognostic_fields(reference_state,
                                                  density,
                                                  momentum,
                                                  potential_temperature_density,
                                                  absolute_humidity,
                                                  condensates,
                                                  tracers)

    timestepper = TimeStepper(timestepper, grid, prognostic_fields)
    pressure_solver = FourierTridiagonalPoissonSolver(grid, tridiagonal_direction=ZDirection())

    model = AtmosphereModel(arch,
                            grid,
                            clock,
                            reference_state,
                            thermodynamics,
                            density,
                            momentum,
                            potential_temperature_density,
                            absolute_humidity,
                            potential_temperature,
                            specific_humidity,
                            temperature,
                            pressure,
                            hydrostatic_pressure_anomaly,
                            pressure_solver,
                            velocities,
                            tracers,
                            advection,
                            coriolis,
                            microphysics,
                            condensates,
                            timestepper)

    update_state!(model)

    return model
end

using Oceananigans.Utils: prettytime, ordered_dict_show, prettykeys
using Oceananigans.TurbulenceClosures: closure_summary

function Base.summary(model::AtmosphereModel)
    A = nameof(typeof(model.grid.architecture))
    G = nameof(typeof(model.grid))
    return string("AtmosphereModel{$A, $G}",
                  "(time = ", prettytime(model.clock.time), ", iteration = ", model.clock.iteration, ")")
end

function Base.show(io::IO, model::AtmosphereModel)
    TS = nameof(typeof(model.timestepper))
    tracernames = prettykeys(model.tracers)

    print(io, summary(model), "\n",
        "├── grid: ", summary(model.grid), "\n",
        "├── timestepper: ", TS, "\n",
        "├── advection scheme: ", summary(model.advection), "\n",
        "├── tracers: ", tracernames, "\n",
        "└── coriolis: ", summary(model.coriolis))
end
