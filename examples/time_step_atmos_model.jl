using AquaSkyLES.AtmosphereModels: AtmosphereModel
using Oceananigans
using Printf

grid = RectilinearGrid(size=(2, 2, 128), x=(0, 1), y=(0, 1), z=(0, 25e3))
model = AtmosphereModel(grid, advection=nothing)

Lz = grid.Lz
Δθ = 5 # K
Tₛ = model.formulation.constants.reference_potential_temperature
θᵢ(x, y, z) = Tₛ + Δθ * z / Lz
qᵢ(x, y, z) = 0
Ξᵢ(x, y, z) = 1e-6 * randn()
set!(model, θ=θᵢ, q=qᵢ, u=Ξᵢ, v=Ξᵢ)

simulation = Simulation(model, Δt=1e-16, stop_iteration=10)

using Printf

ρu, ρv, ρw = model.momentum
δ = Field(∂x(ρu) + ∂y(ρv) + ∂z(ρw))

function progress(sim)
    T = sim.model.temperature
    u, v, w = sim.model.velocities
    T_max = maximum(T)
    T_min = minimum(T)
    u_max = maximum(abs, u)
    v_max = maximum(abs, v)
    w_max = maximum(abs, w)
    compute!(δ)
    δ_max = maximum(abs, δ)
    msg = @sprintf("Iter: %d, time: %s, extrema(T): (%.2f, %.2f) K, max|δ|: %.2e",
                   iteration(sim), prettytime(sim), T_min, T_max, δ_max)
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹", u_max, v_max, w_max)
    @info msg
    return nothing
end

add_callback!(simulation, progress, IterationInterval(1))
