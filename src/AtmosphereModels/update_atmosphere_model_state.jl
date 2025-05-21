#####
##### update state
#####

using Oceananigans.BoundaryConditions: fill_halo_regions!
import Oceananigans.TimeSteppers: update_state!
import Oceananigans: fields, prognostic_fields

const AnelasticModel = AtmosphereModel{<:AnelasticReferenceState}

function prognostic_fields(model::AnelasticModel)
    thermodynamic_fields = (ρθ=model.potential_temperature_density, ρq=model.absolute_humidity)
    return merge(model.momentum, thermodynamic_fields, model.condensates, model.tracers)
end

fields(model::AnelasticModel) = prognostic_fields(model)

function update_state!(model::AnelasticModel)
    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model), async=true)
    compute_auxiliary_variables!(model)
    update_hydrostatic_pressure!(model)
    compute_tendencies!(model)
    return nothing
end

function compute_auxiliary_variables!(model)
    grid = model.grid
    arch = grid.architecture
    velocities = model.velocities
    ref_state = model.reference_state
    momentum = model.momentum

    launch!(arch, grid, :xyz, _compute_velocities!, velocities, ref_state, momentum)

    launch!(arch, grid, :xyz,
            _compute_auxiliary_thermodynamic_variables!,
            model.temperature,
            model.potential_temperature,
            model.specific_humidity,
            model.thermodynamics,
            ref_state,
            model.potential_temperature_density,
            model.absolute_humidity)

    return nothing
end

@kernel function _compute_velocities!(velocities, ref_state, momentum)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        ρʳ = ref_state.density[i, j, k]
        ρu = momentum.ρu[i, j, k]
        ρv = momentum.ρv[i, j, k]
        ρw = momentum.ρw[i, j, k]

        velocities.u[i, j, k] = ρu / ρʳ
        velocities.v[i, j, k] = ρv / ρʳ
        velocities.w[i, j, k] = ρw / ρʳ
    end
end

@kernel function _compute_auxiliary_thermodynamic_variables!(temperature,
                                                             potential_temperature,
                                                             specific_humidity,
                                                             thermo,
                                                             ref_state,
                                                             potential_temperature_density,
                                                             absolute_humidity)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        ρʳ = ref_state.density[i, j, k]
        pʳ = ref_state.pressure[i, j, k]
        ρθ = potential_temperature_density[i, j, k]
        ρq = absolute_humidity[i, j, k]

        potential_temperature[i, j, k] = θ = ρθ / ρʳ
        specific_humidity[i, j, k] = q = ρq / ρʳ
    end

    # Saturation adjustment
    Ψ = AnelasticParcelState(θ, q, ρʳ, pʳ, thermo)
    T = anelastic_temperature(Ψ, thermo)
    @inbounds temperature[i, j, k] = T
end

#####
##### update pressure
#####

using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶜᶜᶠ, ℑzᵃᵃᶠ
using Oceananigans.ImmersedBoundaries: PartialCellBottom, ImmersedBoundaryGrid
using Oceananigans.Grids: topology
using Oceananigans.Grids: XFlatGrid, YFlatGrid

const c = Center()
const f = Face()

@inline function anelastic_buoyancy(i, j, k, grid, ref_state, temperature, specific_humidity, thermo)
    α = anelastic_specific_volume(i, j, k, grid, ref_state, temperature, specific_humidity, thermo)
    αʳ = reference_specific_volume(i, j, k, grid, ref_state, thermo)
    g = thermo.gravitational_acceleration
    return g * (α - αʳ) / αʳ
end

"""
Update the hydrostatic pressure perturbation pHY′. This is done by integrating
the `buoyancy_perturbationᶜᶜᶜ` downwards:

    `pHY′ = ∫ buoyancy_perturbationᶜᶜᶜ dz` from `z=0` down to `z=-Lz`
"""
@kernel function _update_hydrostatic_pressure!(hydrostatic_pressure_anomaly, grid, ref_state, temperature, specific_humidity, thermo)
    i, j = @index(Global, NTuple)
    pₕ′ = hydrostatic_pressure_anomaly

    @inbounds pₕ′[i, j, 1] = anelastic_buoyancy(i, j, 1, grid, ref_state, temperature, specific_humidity, thermo)

    for k in 2:grid.Nz
        Δp′ = ℑzᵃᵃᶠ(i, j, k, grid, anelastic_buoyancy, ref_state, temperature, specific_humidity, thermo)
        @inbounds pₕ′[i, j, k] = pₕ′[i, j, k-1] + Δp′
    end
end

function update_hydrostatic_pressure!(model)
    grid = model.grid
    arch = grid.architecture
    pₕ′ = model.hydrostatic_pressure_anomaly
    ref_state = model.reference_state
    temperature = model.temperature
    specific_humidity = model.specific_humidity
    thermo = model.thermodynamics
    launch!(arch, grid, :xy, _update_hydrostatic_pressure!, pₕ′, grid, ref_state, temperature, specific_humidity, thermo)
    return nothing
end

#####
##### saturation adjustment
#####

function condensate_specific_humidity(T, state::AnelasticParcelState, thermo)
    qᵛ★ = saturation_specific_humidity(T, state.reference_density, thermo, thermo.condensation)
    q = state.specific_humidity
    return max(0, q - qᵛ★)
end

@inline function anelastic_temperature(state::AnelasticParcelState{FT}, thermo) where FT
    θ = state.potential_temperature
    θ == 0 && return zero(FT)

    # Generate guess for unsaturated conditions
    Π = state.exner_function
    T₁ = Π * θ
    qˡ₁ = condensate_specific_humidity(T₁, state, thermo)
    qˡ₁ <= 0 && return T₁
    
    # If we made it this far, we have condensation
    r₁ = saturation_adjustment_residual(T₁, qˡ₁, state, thermo)

    q = state.specific_humidity
    ℒ = thermo.condensation.latent_heat
    cᵖᵐ = mixture_heat_capacity(q, thermo)
    T₂ = (T₁ + sqrt(T₁^2 + 4 * ℒ * qˡ₁ / cᵖᵐ)) / 2
    qˡ₂ = condensate_specific_humidity(T₂, state, thermo)
    r₂ = saturation_adjustment_residual(T₂, qˡ₂, state, thermo)

    # Saturation adjustment
    R = sqrt(max(T₂, T₁))
    ϵ = convert(FT, 1e-4)
    δ = ϵ * R 
    iter = 0

    while abs(r₂ - r₁) > δ
        # Compute slope
        ΔTΔr = (T₂ - T₁) / (r₂ - r₁)

        # Store previous values
        r₁ = r₂
        T₁ = T₂

        # Update
        T₂ -= r₂ * ΔTΔr
        qˡ₂ = condensate_specific_humidity(T₂, state, thermo)
        r₂ = saturation_adjustment_residual(T₂, qˡ₂, state, thermo)
        iter += 1
    end

    return T₂
end

@inline function saturation_adjustment_residual(T, qˡ, state, thermo)
    ℒᵛ = thermo.condensation.latent_heat
    q = state.specific_humidity
    θ = state.potential_temperature
    Π = state.exner_function
    cᵖᵐ = mixture_heat_capacity(q, thermo)
    return T^2 - ℒᵛ * qˡ / cᵖᵐ - Π * θ * T
end

@inline function specific_volume(state, ref, thermo)
    T = temperature(state, ref, thermo)
    Rᵐ = mixture_gas_constant(state.q, thermo)
    pʳ = reference_pressure(state.z, ref, thermo)
    return Rᵐ * T / pʳ
end

function compute_tendencies!(model::AnelasticModel)
    return nothing
end

#=
import Oceananigans.TimeSteppers: update_state!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶜᶜᶠ
using Oceananigans.ImmersedBoundaries: PartialCellBottom, ImmersedBoundaryGrid
using Oceananigans.Grids: topology
using Oceananigans.Grids: XFlatGrid, YFlatGrid

"""
    update_state!(model::AtmosphereModel, callbacks=[])

Update peripheral aspects of the model (halo regions, hydrostatic
pressure) to the current model state. If `callbacks` are provided (in an array),
they are called in the end.
"""
function update_state!(model::AtmosphereModel, callbacks=[]; compute_tendencies = true)

    # Mask immersed tracers
    foreach(model.tracers) do tracer
        mask_immersed_field!(tracer)
    end

    # Fill halos for velocities and tracers
    model_fields = fields(model)
    fill_halo_regions!(model_fields, model.clock, model_fields;
                       fill_boundary_normal_velocities = false, async = true)

    for callback in callbacks
        callback.callsite isa UpdateStateCallsite && callback(model)
    end

    # update_hydrostatic_pressure!(model)
    compute_tendencies && compute_tendencies!(model, callbacks)

    return nothing
end

function compute_tendencies!(model::NonhydrostaticModel, callbacks)

    grid = model.grid
    arch = grid.architecture

    kernel_parameters = :xyz
    compute_interior_tendency_contributions!(model, kernel_parameters; active_cells_map)

    for callback in callbacks
        callback.callsite isa TendencyCallsite && callback(model)
    end

    return nothing
end

""" Store previous value of the source term and compute current source term. """
function compute_interior_tendency_contributions!(model, kernel_parameters; active_cells_map = nothing)

    tendencies           = model.timestepper.Gⁿ
    arch                 = model.architecture
    grid                 = model.grid
    advection            = model.advection
    coriolis             = model.coriolis
    buoyancy             = model.buoyancy
    biogeochemistry      = model.biogeochemistry
    stokes_drift         = model.stokes_drift
    closure              = model.closure
    background_fields    = model.background_fields
    velocities           = model.velocities
    tracers              = model.tracers
    auxiliary_fields     = model.auxiliary_fields
    hydrostatic_pressure = model.pressures.pHY′
    diffusivities        = model.diffusivity_fields
    forcings             = model.forcing
    clock                = model.clock
    u_immersed_bc        = velocities.u.boundary_conditions.immersed
    v_immersed_bc        = velocities.v.boundary_conditions.immersed
    w_immersed_bc        = velocities.w.boundary_conditions.immersed

    start_momentum_kernel_args = (advection,
                                  coriolis,
                                  stokes_drift,
                                  closure)

    end_momentum_kernel_args = (buoyancy,
                                background_fields,
                                velocities,
                                tracers,
                                auxiliary_fields,
                                diffusivities)

    u_kernel_args = tuple(start_momentum_kernel_args...,
                          u_immersed_bc, end_momentum_kernel_args...,
                          hydrostatic_pressure, clock, forcings.u)

    v_kernel_args = tuple(start_momentum_kernel_args...,
                          v_immersed_bc, end_momentum_kernel_args...,
                          hydrostatic_pressure, clock, forcings.v)

    w_kernel_args = tuple(start_momentum_kernel_args...,
                          w_immersed_bc, end_momentum_kernel_args...,
                          hydrostatic_pressure, clock, forcings.w)

    exclude_periphery = true
    launch!(arch, grid, kernel_parameters, compute_Gu!, 
            tendencies.u, grid, u_kernel_args;
            active_cells_map, exclude_periphery)

    launch!(arch, grid, kernel_parameters, compute_Gv!, 
            tendencies.v, grid, v_kernel_args;
            active_cells_map, exclude_periphery)

    launch!(arch, grid, kernel_parameters, compute_Gw!, 
            tendencies.w, grid, w_kernel_args;
            active_cells_map, exclude_periphery)

    start_tracer_kernel_args = (advection, closure)
    end_tracer_kernel_args   = (buoyancy, biogeochemistry, background_fields, velocities,
                                tracers, auxiliary_fields, diffusivities)

    for tracer_index in 1:length(tracers)
        @inbounds c_tendency = tendencies[tracer_index + 3]
        @inbounds forcing = forcings[tracer_index + 3]
        @inbounds c_immersed_bc = tracers[tracer_index].boundary_conditions.immersed
        @inbounds tracer_name = keys(tracers)[tracer_index]

        args = tuple(Val(tracer_index), Val(tracer_name),
                     start_tracer_kernel_args...,
                     c_immersed_bc,
                     end_tracer_kernel_args...,
                     clock, forcing)

        launch!(arch, grid, kernel_parameters, compute_Gc!, 
                c_tendency, grid, args;
                active_cells_map)
    end

    return nothing
end

#####
##### Tendency calculators for u, v, w-velocity
#####

""" Calculate the right-hand-side of the u-velocity equation. """
@kernel function compute_Gu!(Gu, grid, args) 
    i, j, k = @index(Global, NTuple)
    @inbounds Gu[i, j, k] = u_velocity_tendency(i, j, k, grid, args...)
end

""" Calculate the right-hand-side of the v-velocity equation. """
@kernel function compute_Gv!(Gv, grid, args) 
    i, j, k = @index(Global, NTuple)
    @inbounds Gv[i, j, k] = v_velocity_tendency(i, j, k, grid, args...)
end

""" Calculate the right-hand-side of the w-velocity equation. """
@kernel function compute_Gw!(Gw, grid, args) 
    i, j, k = @index(Global, NTuple)
    @inbounds Gw[i, j, k] = w_velocity_tendency(i, j, k, grid, args...)
end


#####
##### Tracer(s)
#####

""" Calculate the right-hand-side of the tracer advection-diffusion equation. """
@kernel function compute_Gc!(Gc, grid, args)
    i, j, k = @index(Global, NTuple)
    @inbounds Gc[i, j, k] = tracer_tendency(i, j, k, grid, args...)
end

#####
##### Boundary contributions to tendencies due to user-prescribed fluxes
#####

""" Apply boundary conditions by adding flux divergences to the right-hand-side. """
function compute_boundary_tendency_contributions!(Gⁿ, arch, velocities, tracers, clock, model_fields)
    fields = merge(velocities, tracers)

    foreach(i -> apply_x_bcs!(Gⁿ[i], fields[i], arch, clock, model_fields), 1:length(fields))
    foreach(i -> apply_y_bcs!(Gⁿ[i], fields[i], arch, clock, model_fields), 1:length(fields))
    foreach(i -> apply_z_bcs!(Gⁿ[i], fields[i], arch, clock, model_fields), 1:length(fields))

    return nothing
end

@kernel function _update_hydrostatic_pressure!(pʰ, grid, buoyancy, C)
    i, j = @index(Global, NTuple)

    @inbounds pʰ[i, j, grid.Nz] = - z_dot_g_bᶜᶜᶠ(i, j, grid.Nz+1, grid, buoyancy, C) * Δzᶜᶜᶠ(i, j, grid.Nz+1, grid)

    for k in grid.Nz-1 : -1 : 1
        @inbounds pʰ[i, j, k] = pʰ[i, j, k+1] - z_dot_g_bᶜᶜᶠ(i, j, k+1, grid, buoyancy, C) * Δzᶜᶜᶠ(i, j, k+1, grid)
    end
end

function update_hydrostatic_pressure!(model; kwargs...)
    grid = model.grid
    arch = model.architecture
    pʰ = model.hydrostatic_pressure_anomaly
    parameters = hydrostatic_pressure_kernel_parameters(grid)
    launch!(arch, grid, parameters, _update_hydrostatic_pressure!, pʰ, grid, model.buoyancy, model.tracers)
end

# extend p kernel to compute also the boundaries
@inline function hydrostatic_pressure_kernel_parameters(grid)
    Nx, Ny, _ = size(grid)
    TX, TY, _ = topology(grid)

    ii = ifelse(TX == Flat, 1:Nx, 0:Nx+1)
    jj = ifelse(TY == Flat, 1:Ny, 0:Ny+1)

    return KernelParameters(ii, jj)
end
=#