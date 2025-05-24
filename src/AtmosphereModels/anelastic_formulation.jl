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

    u = XFaceField(grid)
    v = YFaceField(grid)
    w = ZFaceField(grid)
    velocities = (; u, v, w)

    return velocities, momentum
end

import Oceananigans.Solvers: tridiagonal_direction, _compute_main_diagonal!
import Oceananigans.TimeSteppers: compute_pressure_correction!, make_pressure_correction!

struct AnelasticTridiagonalSolverFormulation{R}
    reference_density :: R
end

tridiagonal_direction(formulation::AnelasticTridiagonalSolverFormulation) = ZDirection()

function formulation_pressure_solver(anelastic_formulation::AnelasticFormulation, grid)
    reference_density = anelastic_formulation.reference_density
    tridiagonal_formulation = AnelasticTridiagonalSolverFormulation(reference_density)
    solver = FourierTridiagonalPoissonSolver(grid; tridiagonal_formulation)
    return solver
end

@kernel function _compute_main_diagonal!(D, grid, λx, λy, formulation::AnelasticTridiagonalSolverFormulation)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)


    # Using a homogeneous Neumann (zero Gradient) boundary condition:
    @inbounds begin
        ρ¹ = @inbounds formulation.reference_density[1, 1, 1]
        ρᴺ = @inbounds formulation.reference_density[1, 1, Nz]

        D[i, j, 1]  = -1 / Δzᵃᵃᶠ(i, j,  2, grid) - ρ¹ * Δzᵃᵃᶜ(i, j,  1, grid) * (λx[i] + λy[j])
        D[i, j, Nz] = -1 / Δzᵃᵃᶠ(i, j, Nz, grid) - ρᴺ * Δzᵃᵃᶜ(i, j, Nz, grid) * (λx[i] + λy[j])

        for k in 2:Nz-1
            ρᵏ = @inbounds formulation.reference_density[1, 1, k]
            ρ⁺ = @inbounds formulation.reference_density[1, 1, k+1]
            ρ⁻ = @inbounds formulation.reference_density[1, 1, k-1]

            ρ̄⁺ = (ρᵏ + ρ⁺) / 2
            ρ̄ᵏ = (ρ⁻ + ρᵏ) / 2

            D[i, j, k] = - (ρ̄⁺ / Δzᵃᵃᶠ(i, j, k+1, grid) + ρ̄ᵏ / Δzᵃᵃᶠ(i, j, k, grid)) - ρᵏ * Δzᵃᵃᶜ(i, j, k, grid) * (λx[i] + λy[j])
        end
    end
end

@kernel function _compute_lower_diagonal!(lower_diagonal, formulation::AnelasticTridiagonalSolverFormulation, grid)
    k = @index(Global)
    @inbounds begin
        ρᵏ = formulation.reference_density[1, 1, k]
        ρ⁺ = formulation.reference_density[1, 1, k+1]
        ρ̄⁺ = (ρᵏ + ρ⁺) / 2
        lower_diagonal[k] = ρ̄⁺ / Δzᵃᵃᶠ(1, 1, k+1, grid)
    end
end

"""
    compute_pressure_correction!(model::NonhydrostaticModel, Δt)

Calculate the (nonhydrostatic) pressure correction associated `tendencies`, `velocities`, and step size `Δt`.
"""
function compute_pressure_correction!(model::AnelasticModel, Δt)
    # Mask immersed velocities
    foreach(mask_immersed_field!, model.momentum)
    fill_halo_regions!(model.momentum, model.clock, fields(model))

    ρʳ = model.formulation.reference_density
    ρU = model.momentum
    solver = model.pressure_solver
    pₙ = model.nonhydrostatic_pressure
    solve_for_anelastic_pressure!(pₙ, solver, ρʳ, ρŨ)

    fill_halo_regions!(pₙ)

    return nothing
end

function solve_for_anelastic_pressure!(pₙ, solver, ρʳ, ρŨ)
    compute_anelastic_source_term!(solver, ρʳ, ρŨ)
    solve!(pₙ, solver)
    return pₙ
end

function compute_anelastic_source_term!(solver::FourierTridiagonalPoissonSolver, ρʳ, ρŨ)
    rhs = solver.source_term
    arch = architecture(solver)
    grid = solver.grid
    launch!(arch, grid, :xyz, _compute_anelastic_source_term!, rhs, grid, ρʳ, ρŨ)
    return nothing
end

@kernel function _compute_anelastic_source_term!(rhs, grid, ρʳ, ρŨ)
    i, j, k = @index(Global, NTuple)
    active = !inactive_cell(i, j, k, grid)
    ρu, ρv, ρw = ρŨ
    δ = divᶜᶜᶜ(i, j, k, grid, ρu, ρv, ρw)
    @inbounds rhs[i, j, k] = active * ρʳ[i, j, k] * Δzᶜᶜᶜ(i, j, k, grid) * δ
end

#=
function compute_source_term!(solver::DistributedFourierTridiagonalPoissonSolver, Ũ)
    rhs = solver.storage.zfield
    arch = architecture(solver)
    grid = solver.local_grid
    tdir = solver.batched_tridiagonal_solver.tridiagonal_direction
    launch!(arch, grid, :xyz, _fourier_tridiagonal_source_term!, rhs, tdir, grid, Ũ)
    return nothing
end
=#

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

function make_pressure_correction!(model::AnelasticModel, Δt)

    launch!(model.architecture, model.grid, :xyz,
            _pressure_correct_momentum!,
            model.momentum,
            model.grid,
            Δt,
            model.nonhydrostatic_pressure)

    return nothing
end