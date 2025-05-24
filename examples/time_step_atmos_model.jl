using AquaSkyLES.AtmosphereModels: AtmosphereModel
using Oceananigans
using Printf

grid = RectilinearGrid(size=(1, 1, 128), x=(0, 1), y=(0, 1), z=(0, 25e3))
model = AtmosphereModel(grid)

Lz = grid.Lz
Δθ = 2 # K
Tₛ = model.formulation.constants.reference_potential_temperature
θᵢ(x, y, z) = Tₛ + Δθ * z / Lz
qᵢ(x, y, z) = 1e-3 + 1e-5 * rand()
Ξᵢ(x, y, z) = 1e-2 * randn()
set!(model, θ=θᵢ, q=qᵢ, u=Ξᵢ, v=Ξᵢ, w=Ξᵢ)

simulation = Simulation(model, Δt=1e-16, stop_iteration=1)

function progress(sim)
    T = sim.model.temperature
    e = sim.model.energy
    q = sim.model.specific_humidity
    ρq = sim.model.absolute_humidity
    u, v, w = sim.model.velocities
    umax, vmax, wmax = maximum(abs, u), maximum(abs, v), maximum(abs, w)
    Tmin, Tmax = minimum(T), maximum(T)
    qmin, qmax = minimum(q), maximum(q)
    msg = @sprintf("Iter: %d, time: %s, extrema(T): (%.2f, %.2f) K",
                   iteration(sim), prettytime(sim), Tmin, Tmax)
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹", umax, vmax, wmax)
    msg *= @sprintf(", extrema(q): (%.2f, %.2f) kg kg⁻¹", qmin, qmax)

    @info msg

    return nothing
end

add_callback!(simulation, progress, IterationInterval(1))

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
