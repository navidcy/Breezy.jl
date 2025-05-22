using ..Thermodynamics:
    saturation_specific_humidity,
    mixture_heat_capacity,
    mixture_gas_constant

using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
import Oceananigans.TimeSteppers: update_state!
import Oceananigans: fields, prognostic_fields

const AnelasticModel = AtmosphereModel{<:AnelasticFormulation}

function prognostic_fields(model::AnelasticModel)
    thermodynamic_fields = (e=model.energy, ρq=model.absolute_humidity)
    return merge(model.momentum, thermodynamic_fields, model.condensates, model.tracers)
end

fields(model::AnelasticModel) = prognostic_fields(model)

function update_state!(model::AnelasticModel, callbacks=[]; compute_tendencies=true)
    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model), async=true)
    compute_auxiliary_variables!(model)
    update_hydrostatic_pressure!(model)
    compute_tendencies && compute_tendencies!(model)
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
            grid,
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
                                                             grid,
                                                             thermo,
                                                             formulation,
                                                             energy,
                                                             absolute_humidity)
    i, j, k = @index(Global, NTuple)

    𝒰 = thermodynamic_state(i, j, k, grid, formulation, thermo, energy, absolute_humidity)
    @inbounds specific_humidity[i, j, k] = 𝒰.specific_humidity

    # Saturation adjustment
    T = compute_temperature(𝒰, thermo)
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

using Oceananigans.Advection: div_𝐯u, div_𝐯v, div_𝐯w, div_Uc
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U, z_f_cross_U
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, ∂zᶜᶜᶠ
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

    scalar_args = (model.advection, model.velocities, model.clock, fields(model))
    Ge = model.timestepper.Gⁿ.e
    e = model.energy
    Fe = model.forcing.e
    e_args = tuple(e, Fe, scalar_args...)
    launch!(arch, grid, :xyz, compute_scalar_tendency!, Ge, grid, e_args)

    ρq = model.absolute_humidity
    Gρq = model.timestepper.Gⁿ.ρq
    Fρq = model.forcing.ρq
    ρq_args = tuple(ρq, Fρq, scalar_args...)
    launch!(arch, grid, :xyz, compute_scalar_tendency!, Gρq, grid, ρq_args)

    return nothing
end

hydrostatic_pressure_gradient_x(i, j, k, grid, pₕ′) = ∂xᶠᶜᶜ(i, j, k, grid, pₕ′)
hydrostatic_pressure_gradient_y(i, j, k, grid, pₕ′) = ∂yᶜᶠᶜ(i, j, k, grid, pₕ′)

@kernel function compute_scalar_tendency!(Gc, grid, args)
    i, j, k = @index(Global, NTuple)
    @inbounds Gc[i, j, k] = scalar_tendency(i, j, k, grid, args...)
end

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

@inline function scalar_tendency(i, j, k, grid,
                                 scalar,
                                 forcing,
                                 advection,
                                 velocities,
                                 clock,
                                 model_fields)

    return ( - div_Uc(i, j, k, grid, advection, velocities, scalar)
             + forcing(i, j, k, grid, clock, model_fields))
end

import Oceananigans.TimeSteppers: calculate_pressure_correction!, pressure_correct_velocities!
using Oceananigans.Models.NonhydrostaticModels: solve_for_pressure!

"""
    calculate_pressure_correction!(model::NonhydrostaticModel, Δt)

Calculate the (nonhydrostatic) pressure correction associated `tendencies`, `velocities`, and step size `Δt`.
"""
function calculate_pressure_correction!(model::AnelasticModel, Δt)

    # Mask immersed velocities
    foreach(mask_immersed_field!, model.momentum)
    fill_halo_regions!(model.momentum, model.clock, fields(model))
    solve_for_pressure!(model.nonhydrostatic_pressure, model.pressure_solver, Δt, model.momentum)
    fill_halo_regions!(model.nonhydrostatic_pressure)

    return nothing
end

#####
##### Fractional and time stepping
#####

"""
Update the predictor velocities u, v, and w with the non-hydrostatic pressure via

    `u^{n+1} = u^n - δₓp_{NH} / Δx * Δt`
"""
@kernel function _pressure_correct_momentum!(M, grid, Δt, pₙ)
    i, j, k = @index(Global, NTuple)

    @inbounds M.ρu[i, j, k] -= ∂xᶠᶜᶜ(i, j, k, grid, pₙ) * Δt
    @inbounds M.ρv[i, j, k] -= ∂yᶜᶠᶜ(i, j, k, grid, pₙ) * Δt
    @inbounds M.ρw[i, j, k] -= ∂zᶜᶜᶠ(i, j, k, grid, pₙ) * Δt
end

function pressure_correct_velocities!(model::AnelasticModel, Δt)

    launch!(model.architecture, model.grid, :xyz,
            _pressure_correct_momentum!,
            model.momentum,
            model.grid,
            Δt,
            model.nonhydrostatic_pressure)

    return nothing
end
