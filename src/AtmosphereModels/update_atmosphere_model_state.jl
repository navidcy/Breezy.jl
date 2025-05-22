using ..MoistThermodynamics:
    saturation_specific_humidity,
    mixture_heat_capacity,
    mixture_gas_constant

using Oceananigans.BoundaryConditions: fill_halo_regions!
import Oceananigans.TimeSteppers: update_state!
import Oceananigans: fields, prognostic_fields

const AnelasticModel = AtmosphereModel{<:AnelasticFormulation}

function prognostic_fields(model::AnelasticModel)
    thermodynamic_fields = (e=model.energy, ρq=model.absolute_humidity)
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
    formulation = model.formulation
    momentum = model.momentum

    launch!(arch, grid, :xyz, _compute_velocities!, velocities, formulation, momentum)

    launch!(arch, grid, :xyz,
            _compute_auxiliary_thermodynamic_variables!,
            model.temperature,
            model.specific_humidity,
            model.thermodynamics,
            formulation,
            model.energy,
            model.absolute_humidity)

    return nothing
end

@kernel function _compute_velocities!(velocities, formulation, momentum)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        ρᵣ = formulation.reference_density[i, j, k]
        ρu = momentum.ρu[i, j, k]
        ρv = momentum.ρv[i, j, k]
        ρw = momentum.ρw[i, j, k]

        velocities.u[i, j, k] = ρu / ρᵣ
        velocities.v[i, j, k] = ρv / ρᵣ
        velocities.w[i, j, k] = ρw / ρᵣ
    end
end

@kernel function _compute_auxiliary_thermodynamic_variables!(temperature,
                                                             specific_humidity,
                                                             thermo,
                                                             formulation,
                                                             energy,
                                                             absolute_humidity)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        ρᵣ = formulation.reference_density[i, j, k]
        pᵣ = formulation.reference_pressure[i, j, k]
        ρq = absolute_humidity[i, j, k]
        e = energy[i, j, k]

        specific_humidity[i, j, k] = q = ρq / ρᵣ
    end

    # Saturation adjustment
    cₚ = thermo.dry_air.heat_capacity
    θ = e / (cₚ * ρᵣ)
    p₀ = formulation.constants.base_pressure
    Ψ = AnelasticThermodynamicState(θ, q, ρᵣ, pᵣ, p₀, thermo)
    T = compute_temperature(Ψ, thermo)
    @inbounds temperature[i, j, k] = T
end

#=
@inline function specific_volume(state, ref, thermo)
    T = temperature(state, ref, thermo)
    Rᵐ = mixture_gas_constant(state.q, thermo)
    pᵣ = reference_pressure(state.z, ref, thermo)
    return Rᵐ * T / pᵣ
end
=#

using Oceananigans.Utils: launch!

function compute_tendencies!(model::AnelasticModel)
    grid = model.grid
    arch = grid.architecture
    Gρu = model.timestepper.Gⁿ.ρu
    Gρv = model.timestepper.Gⁿ.ρv
    Gρw = model.timestepper.Gⁿ.ρw

    common_args = (model.advection,
                   model.velocities,
                   model.momentum,
                   model.coriolis,
                   model.clock,
                   fields(model))    

    pₕ′ = model.hydrostatic_pressure_anomaly
    u_args = tuple(common_args..., model.forcing.ρu, pₕ′)
    v_args = tuple(common_args..., model.forcing.ρv, pₕ′)
    w_args = tuple(common_args..., model.forcing.ρw)

    launch!(arch, grid, :xyz, compute_x_momentum_tendency!, Gρu, grid, u_args)
    launch!(arch, grid, :xyz, compute_y_momentum_tendency!, Gρv, grid, v_args)
    launch!(arch, grid, :xyz, compute_z_momentum_tendency!, Gρw, grid, w_args)

    return nothing
end

using Oceananigans.Advection: div_𝐯u, div_𝐯v, div_𝐯w
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U, z_f_cross_U
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ

hydrostatic_pressure_gradient_x(i, j, k, grid, pₕ′) = ∂xᶠᶜᶜ(i, j, k, grid, pₕ′)
hydrostatic_pressure_gradient_y(i, j, k, grid, pₕ′) = ∂yᶜᶠᶜ(i, j, k, grid, pₕ′)

@kernel function compute_x_momentum_tendency!(Gρu, grid, args)
    i, j, k = @index(Global, NTuple)
    @inbounds Gρu[i, j, k] = x_momentum_tendency(i, j, k, grid, args...)
end

@kernel function compute_y_momentum_tendency!(Gρv, grid, args)
    i, j, k = @index(Global, NTuple)
    @inbounds Gρv[i, j, k] = y_momentum_tendency(i, j, k, grid, args...)
end

@kernel function compute_z_momentum_tendency!(Gρw, grid, args)
    i, j, k = @index(Global, NTuple)
    @inbounds Gρw[i, j, k] = z_momentum_tendency(i, j, k, grid, args...)
end

@inline function x_momentum_tendency(i, j, k, grid,
                                     advection,
                                     velocities,
                                     momentum,
                                     coriolis,
                                     clock,
                                     model_fields,
                                     forcing,
                                     hydrostatic_pressure_anomaly)

    return ( - div_𝐯u(i, j, k, grid, advection, velocities, momentum.ρu)
             - x_f_cross_U(i, j, k, grid, coriolis, momentum)
             - hydrostatic_pressure_gradient_x(i, j, k, grid, hydrostatic_pressure_anomaly)
             + forcing(i, j, k, grid, clock, model_fields))
end

@inline function y_momentum_tendency(i, j, k, grid,
                                     advection,
                                     velocities,
                                     momentum,
                                     coriolis,
                                     clock,
                                     model_fields,
                                     forcing,
                                     hydrostatic_pressure_anomaly)

    return ( - div_𝐯v(i, j, k, grid, advection, velocities, momentum.ρu)
             - y_f_cross_U(i, j, k, grid, coriolis, momentum)
             - hydrostatic_pressure_gradient_y(i, j, k, grid, hydrostatic_pressure_anomaly)
             + forcing(i, j, k, grid, clock, model_fields))
end

@inline function z_momentum_tendency(i, j, k, grid,
                                     advection,
                                     velocities,
                                     momentum,
                                     coriolis,
                                     clock,
                                     model_fields,
                                     forcing)

    return ( - div_𝐯v(i, j, k, grid, advection, velocities, momentum.ρu)
             - z_f_cross_U(i, j, k, grid, coriolis, momentum)
             + forcing(i, j, k, grid, clock, model_fields))
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