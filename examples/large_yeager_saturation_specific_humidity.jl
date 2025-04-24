using GLMakie
using AquaSkyLES

thermo = AquaSkyLES.AtmosphereThermodynamics()

saturation_specific_humidity_large_yeager(T, ρ) = 640380 * exp(-5107.4 / T) / ρ

T₀ = 273.15
T = collect(T₀:0.01:T₀+50)
q★_large_yeager = saturation_specific_humidity_large_yeager.(T, 1.2)
q★_aqua_sky = AquaSkyLES.saturation_specific_humidity.(T, 1.2, Ref(thermo))

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Temperature (K)", ylabel = "Saturation Specific Humidity (g/kg)")
lines!(ax, T, q★_large_yeager, label = "Large and Yeager (2009)")
lines!(ax, T, q★_aqua_sky, label = "AquaSkyLES")
Legend(fig[1, 2], ax)
display(fig)