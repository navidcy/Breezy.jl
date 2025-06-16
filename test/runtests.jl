using Test
using Breeze

@testset "Breeze.jl" begin
    @testset "Thermodynamics" begin
        thermo = Breeze.AtmosphereThermodynamics()

        # Test saturation specific humidity calculation
        T = 293.15  # 20°C
        ρ = 1.2     # kg/m³
        q★ = Breeze.saturation_specific_humidity(T, ρ, thermo, thermo.condensation)

        # Basic sanity checks
        @test q★ > 0
        @test q★ < 1  # specific humidity should be less than 1

        # Test temperature dependence
        T_cold = 273.15  # 0°C
        q★_cold = Breeze.saturation_specific_humidity(T_cold, ρ, thermo, thermo.condensation)
        @test q★ > q★_cold  # warmer air should hold more moisture
    end
end
