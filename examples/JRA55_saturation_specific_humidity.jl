using JLD2
using AquaSkyLES: density, saturation_specific_humidity, MoistAuxiliaryState, AtmosphereThermodynamics
using GLMakie

@load "JRA55_atmospheric_state_Jan_1_1991.jld2" q T p

thermo = AtmosphereThermodynamics()

A = @. MoistAuxiliaryState(p, T, q)
ρ = density.(A, Ref(thermo))
q★ = saturation_specific_humidity.(A, Ref(thermo))

scatter(q★, T, color=ρ, colormap=:blues)