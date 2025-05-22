using AquaSkyLES.AtmosphereModels: AtmosphereModel
using Oceananigans

grid = RectilinearGrid(size=(1, 1, 128), x=(0, 1), y=(0, 1), z=(0, 25e3))
model = AtmosphereModel(grid)

cᵖ = model.thermodynamics.dry_air.heat_capacity
ρ₀ = first(model.formulation.reference_density)
Lz = grid.Lz
Δθ = 5 # K
Tₛ = model.formulation.constants.reference_potential_temperature
eᵢ(x, y, z) = ρ₀ * cᵖ * (Tₛ + Δθ * z / Lz)
qᵢ(x, y, z) = 1e-2 + 1e-5 * rand()
set!(model, e=eᵢ, q=qᵢ)

using GLMakie

using AquaSkyLES.MoistThermodynamics: dry_air_gas_constant
Rᵈ = dry_air_gas_constant(model.thermodynamics)
ϰ = cᵖ / Rᵈ
pᵣ = model.formulation.reference_pressure
Π = Field((pᵣ / pᵣ[1])^(1/ϰ))
compute!(Π)

θ = CenterField(grid)
set!(θ, (x, y, z) -> Tₛ + Δθ * z / Lz)
T1 = Field(θ * Π)
compute!(T1)

axT = Axis(fig[1, 1], title="Temperature")
#lines!(axT, model.temperature)
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
