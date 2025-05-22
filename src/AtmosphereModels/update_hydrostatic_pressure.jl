#####
##### update pressure
#####

using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶜᶜᶠ, ℑzᵃᵃᶠ
using Oceananigans.ImmersedBoundaries: PartialCellBottom, ImmersedBoundaryGrid
using Oceananigans.Grids: topology
using Oceananigans.Grids: XFlatGrid, YFlatGrid

const c = Center()
const f = Face()

@inline function buoyancy(i, j, k, grid, formulation, temperature, specific_humidity, thermo)
    α = specific_volume(i, j, k, grid, formulation, temperature, specific_humidity, thermo)
    αʳ = reference_specific_volume(i, j, k, grid, formulation, thermo)
    g = thermo.gravitational_acceleration
    return g * (α - αʳ) / αʳ
end

"""
Update the hydrostatic pressure perturbation pHY′. This is done by integrating
the `buoyancy_perturbationᶜᶜᶜ` downwards:

    `pHY′ = ∫ buoyancy_perturbationᶜᶜᶜ dz` from `z=0` down to `z=-Lz`
"""
@kernel function _update_hydrostatic_pressure!(pₕ′, grid, formulation, ρᵣ, T, q, thermo)
    i, j = @index(Global, NTuple)

    b₁ = buoyancy(i, j, 1, grid, formulation, T, q, thermo)
    @inbounds pₕ′[i, j, 1] = ρᵣ[i, j, 1] * b₁

    @inbounds for k in 2:grid.Nz
        bₖ = ℑzᵃᵃᶠ(i, j, k, grid, buoyancy, formulation, T, q, thermo)
        Δp′ = ρᵣ[i, j, k] * bₖ
        pₕ′[i, j, k] = pₕ′[i, j, k-1] + Δp′
    end
end

function update_hydrostatic_pressure!(model)
    grid = model.grid
    arch = grid.architecture
    pₕ′ = model.hydrostatic_pressure_anomaly
    formulation = model.formulation
    T = model.temperature
    q = model.specific_humidity
    ρᵣ = model.formulation.reference_density
    thermo = model.thermodynamics
    launch!(arch, grid, :xy, _update_hydrostatic_pressure!, pₕ′, grid, formulation, ρᵣ, T, q, thermo)
    return nothing
end
