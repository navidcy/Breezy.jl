using Oceananigans
using Oceananigans.Models: AbstractModel
using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.TimeSteppers: TimeStepper

struct BoussinesqReferenceState{FT}
    base_density :: FT
end

field_names(::BoussinesqReferenceState, tracer_names) = (:u, :v, :w, tracer_names...)

function make_prognostic_fields(::BoussinesqReferenceState, density, velocities, momentum, total_energy, specific_humidity, condensates, tracers)
    merge(velocities,
          (e=total_energy, q=specific_humidity),
          condensates,
          tracers)
end

function default_reference_state(grid)
    FT = eltype(grid)
    return BoussinesqReferenceState{FT}(1.2)
end

function materialize_momentum_and_velocities(reference_state::BoussinesqReferenceState, grid, boundary_conditions)
    u = XFaceField(grid, boundary_conditions=boundary_conditions.u)
    v = YFaceField(grid, boundary_conditions=boundary_conditions.v)
    w = ZFaceField(grid, boundary_conditions=boundary_conditions.w)
    velocities = (; u, v, w)

    ρ₀ = reference_state.base_density
    ρu = ρ₀ * u
    ρv = ρ₀ * v
    ρw = ρ₀ * w
    momentum = (; ρu, ρv, ρw)

    return velocities, momentum
end

struct WarmPhaseSaturationAdjustment end

struct AtmosphereModel{Ts, Re, Ar, Gr, Cl, Ad, De, Mo, En, Hu, Te, Pr, Pp, Ve, Tr, Co, Mi, Cn} <: AbstractModel{Ts, Ar}
    architecture :: Ar
    grid :: Gr
    clock :: Cl
    reference_state :: Re
    density :: De
    momentum :: Mo
    total_energy :: En
    specific_humidity :: Hu
    temperature :: Te
    pressure :: Pr
    hydrostatic_pressure_anomaly :: Pp
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
                         reference_state = default_reference_state(grid),
                         total_energy = CenterField(grid),
                         specific_humidity = CenterField(grid),
                         coriolis = nothing,
                         boundary_conditions = NamedTuple(),
                         advection = WENO(order=5),
                         microphysics = WarmPhaseSaturationAdjustment(),
                         timestepper = :RungeKutta3)

    arch = architecture(grid)
    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)

    hydrostatic_pressure_anomaly = CenterField(grid)
    pressure = CenterField(grid)

    # Next, we form a list of default boundary conditions:
    names = field_names(reference_state, tracers)
    default_boundary_conditions = NamedTuple{names}(FieldBoundaryConditions() for _ in names)
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, names)

    velocities, momentum = materialize_momentum_and_velocities(reference_state, grid, boundary_conditions)
    tracers = NamedTuple(n => CenterField(grid, boundary_conditions=boundary_conditions[n]) for name in tracers)
    condensates = materialize_condenstates(microphysics, grid)
    advection = adapt_advection_order(advection, grid)

    prognostic_fields = make_prognostic_fields(reference_state, density, velocities, momentum, total_energy, tracers)
    timestepper = TimeStepper(timestepper, grid, prognostic_fields)
    
    model = AtmosphereModel(arch,
                            grid,
                            clock,
                            reference_state,
                            density,
                            momentum,
                            total_energy,
                            specific_humidity,
                            temperature,
                            pressure,
                            velocities,
                            tracers,
                            advection,
                            coriolis,
                            microphysics,
                            condensates,
                            timestepper)

    return model
end
   
