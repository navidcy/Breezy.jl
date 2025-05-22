using AquaSkyLES.AtmosphereModels: AtmosphereModel
using Oceananigans

grid = RectilinearGrid(size=(1, 1, 128), x=(0, 1), y=(0, 1), z=(0, 25e3))

model = AtmosphereModel(grid)

Lz = grid.Lz
Δθ = 1 # K
Tₛ = model.formulation.constants.reference_potential_temperature
θᵢ(x, y, z) = Tₛ + Δθ * z / Lz
qᵢ(x, y, z) = 1e-2
Ξᵢ(x, y, z) = 1e-2 * randn()
set!(model, θ=θᵢ, q=qᵢ, u=Ξᵢ, v=Ξᵢ, w=Ξᵢ)

simulation = Simulation(model, Δt=10, stop_iteration=1)
run!(simulation)

using GLMakie

using AquaSkyLES.Thermodynamics: dry_air_gas_constant
cᵖᵈ = model.thermodynamics.dry_air.heat_capacity
Rᵈ = dry_air_gas_constant(model.thermodynamics)
pᵣ = model.formulation.reference_pressure
Π = Field((pᵣ / pᵣ[1, 1, 1])^(Rᵈ / cᵖᵈ))
compute!(Π)

θ = CenterField(grid)
set!(θ, (x, y, z) -> Tₛ + Δθ * z / Lz)
T1 = Field(θ * Π)
compute!(T1)

fig = Figure(size=(1000, 1000))

axT = Axis(fig[1, 1], title="Temperature")
lines!(axT, model.temperature)
lines!(axT, T1)

axp = Axis(fig[1, 2], title="Reference Pressure")
lines!(axp, model.formulation.reference_pressure)

axp = Axis(fig[1, 3], title="Reference Density")
lines!(axp, model.formulation.reference_density)

axq = Axis(fig[1, 4], title="Specific Humidity")
lines!(axq, model.specific_humidity)

axe = Axis(fig[1, 5], title="Energy")
lines!(axe, model.energy)

fig
